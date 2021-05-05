using Random: AbstractRNG, Random.GLOBAL_RNG

using DataStructures: BinaryMinHeap
using Distributions: Binomial, Discrete, Gamma, Sampleable, Univariate
using FileIO: load
using Images: Gray
using NamedDims: NamedDimsArray
using YAML: load_file

const HOURS_IN_DAY = 24
const DAYS_IN_WEEK = 7
const TIME_HORIZON = 840
hour(time::Integer)::Int = (time - 1) % HOURS_IN_DAY + 1
day(time::Integer)::Int = ((time - 1) ÷ HOURS_IN_DAY) % DAYS_IN_WEEK + 1
is_weekend(time::Integer)::Bool = day(time) == 6 || day(time) == 7

"""
    Strain(initial, strength, radius, incubation_mean, incubation_shape, infection_mean, infection_shape)

Stores the parameters describing the dynamics of a virus strain.

**Fields*
- `initial (Float)`: The probability of initial infection.
- `strength (Float64)`: The infection strength of the virus strain.
- `radius (Float64)`: The maximum infection radius of the virus strain (in metres).
- `incubation_mean (Float64)`: The mean parameter used by the gamma-distributed incubation
    durations.
- `incubation_shape (Float64)`: The shape parameter used by the gamma-distributed incubation
    durations.
- `infection_mean (Float64)`: The mean parameter used by the gamma-distributed infection
    durations.
- `infection_shape (Float64)`: The shape parameter used by the gamma-distributed infection
    durations.
"""
struct Strain
    initial::Float64
    strength::Float64
    radius::Float64
    incubation_mean::Float64
    incubation_shape::Float64
    infection_mean::Float64
    infection_shape::Float64
end

"""
    Strains(initial, strength, radius, infection_mean, infection_shape)

Stores the parameters describing the dynamics of multiple virus strains.

**Fields*
- `initial (Vector{Float64}): The probabilities of initial infection.
- `strength (Vector{Float64})`: The infection strengths of the virus strains.
- `radius (Vector{Float64})`: The maximum infection radii of the virus strains (in metres).
- `infection_mean (Vector{Float64})`: The mean infection duration.
- `infection_shape (Vector{Float64})`: The shape parameters used by the gamma-distributed
    infection durations.

Alternatively, when any argument is replaced by a `Float64` instead of `Vector{Float64}`,
the underlying vector will be created automatically.
"""
struct Strains
    initial::Vector{Float64}
    strength::Vector{Float64}
    radius::Vector{Float64}
    infection_mean::Vector{Float64}
    infection_shape::Vector{Float64}
end

function Strains(
    initial::Union{Real, Vector{T}},
    strength::Union{Real, Vector{T}},
    radius::Union{Real, Vector{T}},
    infection_mean::Union{Real, Vector{T}},
    infection_shape::Union{Real, Vector{T}}
) where T <: Real
    make_vector(x) = Float64.(x isa Real ? [x] : x)

    return Strains(
        make_vector(initial),
        make_vector(strength),
        make_vector(radius),
        make_vector(infection_mean),
        make_vector(infection_shape)
    )
end

"""
    State(time, susceptible, exposed, infected, recovered, infection_times, recovery_times)

Represents the epidemic state of a virus strain within a population.

**Fields**
- `time::Int`: The current time (measured in hours from 12:00am Monday).
- `susceptible::Int`: The number of susceptible individuals.
- `exposed::Int`: The number of exposed individuals.
- `infected::Int`: The number of infected individuals.
- `recovered::Int`: The number of recovered individuals.
- `infection_times::Vector{Float64}`: A collection of infection times.
- `recovery_times::Vector{Float64}`: A collection of recovery times.


    State(rng, population, strain)

Generates the initial state of a vrius strain within a population.

**Fields**
- `rng::AbstractRNG`: A random number generator.
- `population::Integer`: The initial number of participants.
- `strain::Strain`: Stores the parameters describing the virus strain.
"""
mutable struct State
    time::Int
    susceptible::Int
    exposed::Int
    infected::Int
    recovered::Int
    infection_times::BinaryMinHeap{Float64}
    recovery_times::BinaryMinHeap{Float64}
end

function State(rng::AbstractRNG, population::Integer, strain::Strain)
    infected = rand(Binomial(population, strain.initial))
    state = State(
        0, population - infected, 0, infected, 0, BinaryMinHeap{Float64}(),
        BinaryMinHeap{Float64}()
    )

    for _ in 1:infected
        schedule_recovery!(rng, state, strain)
    end

    return state
end

"""
    CampusSampler(rows, columns, scale, probabilities, aliases)

A sampler used for generating individuals' locations using the alias method.

**Fields**
- `rows::Int`: The number of rows used to represent the campus.
- `columns::Int`: The number of columns used to represent the campus.
- `scale::Float64`: The scale used when representing the campus as a matrix
    (in metres/cell).
- `probabilities::Matrix{Float64}`: The probability table used by the alias method when
    sampling.
- `aliases::Matrix{Int}`: The alias table used by the alias method when sampling.


    CampusSampler(scale, weights)

Builds a `CampusSampler` from a matrix of weights.
"""
struct CampusSampler <: Sampleable{Univariate, Discrete}
    rows::Int
    columns::Int
    scale::Float64
    probabilities::Matrix{Float64}
    aliases::Matrix{Int}
end

function CampusSampler(scale::Float64, weights::Matrix{Float64})
    rows, columns = size(weights)

    probabilities = weights / sum(weights) * rows * columns
    aliases = zeros(Int, size(weights))

    indices = LinearIndices(size(weights))
    under = findall(x -> x < 1.0, probabilities)
    over = findall(x -> x >= 1.0, probabilities)

    while length(over) != 0 && length(under) != 0
        α, β = pop!(under), over[end]
        aliases[α] = indices[β]
        probabilities[β] += probabilities[α] - 1

        if probabilities[β] < 1.0
            pop!(over)
            push!(under, β)
        end
    end

    return CampusSampler(rows, columns, scale, probabilities, aliases)
end

"""
    Base.rand(rng, sampler)

Returns a randomly generated position of an individual on campus.

**Arguments**
- `rng::AbstractRNG`: A random number generator.
- `sampler::CampusSampler`: A sampler describing the distribution of positions on campus.
"""
function Base.rand(rng::AbstractRNG, sampler::CampusSampler)
    u = sampler.rows * sampler.columns * rand(rng)
    k = Int(floor(u)) + 1

    if u - k + 1 >= sampler.probabilities[k]
        k = sampler.aliases[k]
    end

    x = sampler.scale * (((k - 1) % sampler.rows) + rand(rng))
    y = sampler.scale * (((k - 1) ÷ sampler.rows) + rand(rng))

    return (x=x, y=y)
end

"""
    Behaviour(campus_sampler, weekday_attendance, weekend_attendance, compliance)

Stores the population dynamics of a campus.

**Fields**
- `campus_sampler::CampusSampler`: Generates the positions of individuals on campus.
- `weekday_attendance::Vector{Float64}`: The attendance probabilities for each hour of a
    weekday.
- `weekend_attendance::Vector{Float64}`: The attendance probabilities for each hour of a
    weekend.
- `compliance::Float64`: The probability that an individual is running the Safe Blues app.
"""
struct Behaviour
    campus_sampler::CampusSampler
    weekday_attendance::Vector{Float64}
    weekend_attendance::Vector{Float64}
    compliance::Float64
end

"""
    load_behaviour(behaviour_file, heatmap_file)

Loads an `Behaviour` from an behaviour configuration and a campus heatmap.

**Arguments**
- `behaviour_file::String`: The path of an behaviour configuration file (`.yaml`).
- `heatmap_file::String`: The path of a campus heatmap file (`.png`).
"""
function load_behaviour(behaviour_file::String, heatmap_file::String)
    contents = load_file(behaviour_file)

    scale::Float64 = contents["scale"]
    weekday_attendance::Vector{Float64} = contents["weekday_attendance"]
    weekend_attendance::Vector{Float64} = contents["weekend_attendance"]
    compliance::Float64 = contents["compliance"]

    weights = (x -> Float64(Gray(x))).(load(heatmap_file))
    campus_sampler = CampusSampler(scale, weights)

    return Behaviour(campus_sampler, weekday_attendance, weekend_attendance, compliance)
end

const GLOBAL_BEHAVIOUR = cd(@__DIR__) do
    return load_behaviour("assets/behaviour.yaml", "assets/heatmap.png")
end

"""
    Intervention(start, stop, strength)

Stores the parameters describing a social distancing intervention.

**Fields**
- `start::Int`: The first time period (in hours from 12:00am Monday) of the intervention.
- `stop::Int`: The last time period (in hours from 12:00am Monday) of the intervention.
- `strength::Float64`: The strength of the intervention (from weakest `0.0` to strongest
    `1.0`).
"""
struct Intervention
    start::Int
    stop::Int
    strength::Float64
end

const DEFAULT_INTERVENTION = Intervention(0, 0, 0.0)

"""
    schedule_infection!(rng, state, strain)

Generates the infection time of a newly exposed individual.

Inserts the new infection time into `state.infection_times`.

**Arguments**
- `rng::AbstractRNG`: A random number generator.
- `state::State`: The epidemic state of the strain.
- `strain::Strain`: Stores the parameters describing the virus strain.
"""
function schedule_infection!(rng::AbstractRNG, state::State, strain::Strain)
    time = state.time + rand(rng, Gamma(
        strain.incubation_shape, strain.incubation_mean / strain.incubation_shape
    ))

    # Ignore infinite incubation durations. This individual will never become infected.
    if time == Inf
        return state
    end

    # Insert the infection time.
    push!(state.infection_times, time)

    return state
end

"""
    schedule_recovery!(rng, state, strain)

Generates the recovery time of a newly infected individual.

Inserts the new recovery time into `state.recovery_times`.

**Arguments**
- `rng::AbstractRNG`: A random number generator.
- `state::State`: The epidemic state of the strain.
- `strain::Strain`: Stores the parameters describing the virus strain.
"""
function schedule_recovery!(rng::AbstractRNG, state::State, strain::Strain)
    time = state.time + rand(rng, Gamma(
        strain.infection_shape, strain.infection_mean / strain.infection_shape
    ))

    # Ignore infinite infections durations. This individual will never recover.
    if time == Inf
        return state
    end

    # Insert the recovery time.
    push!(state.recovery_times, time)

    return state
end

"""
    spread(rng, state, strain, behaviour, intervention)

Returns the number of newly exposed individuals in a single time period.

**Arguments**
- `rng::AbstractRNG`: A random number generator.
- `state::State`: The epidemic state of the strain.
- `strain::Strain`: Stores the parameters describing the virus strain.
- `behaviour::Behaviour`: Stores the population dynamics of the campus.
- `intervention::Intervention`: Stores the parameters describing a social distancing
    intervention.
"""
function spread(
    rng::AbstractRNG,
    state::State,
    strain::Strain,
    behaviour::Behaviour,
    intervention::Intervention
)
    # Get the number of active (attending & complying) participants.
    activity = (
        is_weekend(state.time) ? behaviour.weekend_attendance
        : behaviour.weekday_attendance
    )[hour(state.time)] * behaviour.compliance
    active_susceptible = rand(rng, Binomial(state.susceptible, activity))
    active_infected = rand(rng, Binomial(state.infected, activity))
    if active_susceptible == 0 || active_infected == 0
        return 0
    end

    # Calculate the adjusted infection strength.
    intervene = intervention.start <= state.time <= intervention.stop
    strength = strain.strength * (intervene ? 1.0 - intervention.strength : 1.0)

    # Generate the locations of infected individuals.
    infected_points = [rand(rng, behaviour.campus_sampler) for _ in 1:active_infected]

    bound = strain.radius^2

    # Attempt infections between nearby individuals.
    exposed = 0
    for _ in 1:active_susceptible
        point = rand(rng, behaviour.campus_sampler)

        p = 0
        for point′ in infected_points
            distance = (point.x - point′.x)^2 + (point.y - point′.y)^2
            if distance >= bound
                continue
            end

            q = 1 - exp(-strength * (1 - √distance / strain.radius))
            p = p + q - p * q
        end
    
        # Determine whether an infection occurs.
        if rand(rng) < p
            exposed += 1
        end
    end

    return exposed
end

"""
    advance!(rng, state, strain, behaviour, intervention)

Advances the epidemic state of a virus strain forward by a single time increment.

Returns a `NamedTuple` storing the number of susceptible (`.susceptible`), infected
(`.infected`), and recovered (`.recovered`) individuals. Updates `state` to reflect the new
epidemic state of the virus.

**Arguments**
- `rng::AbstractRNG`: A random number generator.
- `state::State`: The epidemic state of the strain.
- `strain::Strain`: Stores the parameters describing the virus strain.
- `behaviour::Behaviour`: Stores the population dynamics of the campus.
- `intervention::Intervention`: Stores the parameters describing a social distancing
    intervention.

**Keyword Arguments**
- `mode::Symbol=:SIR`: Selects the epidemic model used from one of :SIR, :SIS, :SI, :SEIR,
    :SEIS, or :SEI.
"""
function advance!(
    rng::AbstractRNG,
    state::State,
    strain::Strain,
    behaviour::Behaviour,
    intervention::Intervention;
    mode::Symbol=:SIR
)
    state.time += 1
    susceptible_in, exposed_in, infected_in, recovered_in = 0, 0, 0, 0

    # Handle movement out of susceptible compartment.
    susceptible_out = spread(rng, state, strain, behaviour, intervention)
    state.susceptible -= susceptible_out

    if mode == :SIR || mode == :SIS || mode == :SI
        infected_in += susceptible_out
    elseif mode == :SEIR || mode == :SEIS || mode == :SEI
        exposed_in += susceptible_out
    else
        throw(ArgumentError("invalid mode '$mode'"))
    end

    # Handle movement into exposed compartment.
    state.exposed += exposed_in
    for _ in 1:exposed_in
        schedule_infection!(rng, state, strain)
    end

    # Handle movement out of exposed compartment.
    exposed_out = 0
    while length(state.infection_times) != 0
        time = first(state.infection_times)
        if time <= state.time
            pop!(state.infection_times)
            exposed_out += 1
            continue
        end

        break
    end

    infected_in += exposed_out
    state.exposed -= exposed_out

    # Handle movement into infected compartment.
    state.infected += infected_in
    for _ in 1:infected_in
        schedule_recovery!(rng, state, strain)
    end

    # Handle movement out of infected compartment.
    infected_out = 0
    while length(state.recovery_times) != 0
        time = first(state.recovery_times)
        if time <= state.time
            pop!(state.recovery_times)
            infected_out += 1
            continue
        end

        break
    end

    state.infected -= infected_out

    if mode == :SIR || mode == :SEIR
        state.recovered += infected_out
    elseif mode == :SIS || mode == :SEIS
        state.susceptible += infected_out
    elseif mode == :SI || mode == :SEI
        state.infected += infected_out
    else
        throw(ArgumentError("invalid mode '$mode'"))
    end

    return (
        susceptible=state.susceptible,
        exposed=state.exposed,
        infected=state.infected,
        recovered=state.recovered
    )
end

SimulationArray = NamedDimsArray{(:time, :trial)}
SimulationData = @NamedTuple begin
    population::Int
    susceptible::SimulationArray
    exposed::SimulationArray
    infected::SimulationArray
    recovered::SimulationArray
end

ParametricArray = NamedDimsArray{(:time, :trial, :initial, :strength, :radius, :infection_mean, :infection_shape)}
ParametricData = @NamedTuple begin
    population::Int
    strains::Strains

    susceptible::ParametricArray
    infected::ParametricArray
    recovered::ParametricArray
end

"""
    simulate(rngs, population, strain, behaviour)

Simulates the spread of a virus strain within a population and returns `SimulationData`.

**Arguments**
- `rngs::Vector{AbstractRNG}`: A random number generator for each simulation trial.
- `population::Integer`: The initial number of participants.
- `strain::Strain`: Stores the parameters describing the virus strain.
- `behaviour::Behaviour`: Stores the population dynamics of the campus.

**Keyword Arguments**
- `intervention::Intervention=DEFAULT_INTERVENTION`: Stores the parameters describing a
    social distancing intervention.
- `mode::Symbol=:SIR`: Selects the epidemic model used from one of :SIR, :SIS, :SI, :SEIR,
    :SEIS, or :SEI.
- `arrivals::Integer=0`: The number of new individuals that will arrive throughout the
    simulation period.


    simulate(population, strain, behaviour)
    simulate(rngs, population, strain)
    simulate(population, strain)

Alternatively, when `rngs` is ommited the random number generators default to
`[Random.GLOBAL_RNG]` and when `behaviour` is ommited the behaviour defaults to
`GLOBAL_BEHAVIOUR`.


    simulate(rng, population, strain, behaviour)
    simulate(rng, population, strain)

If only a single random number generator `rng::AbstractRNG` is passed to `simulate`, then
the simulation will be run exactly once.


    simulate(rng, population, strains, behaviour)

Finally, if `strain::Strain` is replaced by `strains::Strains`, then multiple simulations
will be run using different strain parameters and `ParametricData` is returned.
"""
function simulate(
    rngs::Vector{T},
    population::Integer,
    strain::Strain,
    behaviour::Behaviour;
    intervention::Intervention=DEFAULT_INTERVENTION,
    mode::Symbol=:SIR,
    arrivals::Integer=0
)::SimulationData where T <: AbstractRNG
    trials = length(rngs)
    susceptible = SimulationArray(zeros(Int, TIME_HORIZON + 1, trials))
    exposed = SimulationArray(zeros(Int, TIME_HORIZON + 1, trials))
    infected = SimulationArray(zeros(Int, TIME_HORIZON + 1, trials))
    recovered = SimulationArray(zeros(Int, TIME_HORIZON + 1, trials))

    arrival_rate = arrivals / TIME_HORIZON
    arrival_excess = 0.0

    rngs = copy.(rngs)

    for trial in 1:trials
        state = State(rngs[trial], population, strain)

        susceptible[time=1, trial=trial] = state.susceptible
        exposed[time=1, trial=trial] = state.exposed
        infected[time=1, trial=trial] = state.infected
        recovered[time=1, trial=trial] = state.recovered

        for time in 2:(TIME_HORIZON + 1)
            # Calculate the number of new arrivals.
            arrival_excess += arrival_rate
            state.susceptible += floor(arrival_excess)
            arrival_excess %= 1.0

            advance!(rngs[trial], state, strain, behaviour, intervention; mode=mode)

            susceptible[time=time, trial=trial] = state.susceptible
            exposed[time=time, trial=trial] = state.exposed
            infected[time=time, trial=trial] = state.infected
            recovered[time=time, trial=trial] = state.recovered
        end
    end

    return (
        population=population + arrivals,
        susceptible=susceptible,
        exposed=exposed,
        infected=infected,
        recovered=recovered
    )
end

function simulate(
    rngs::Vector{T},
    population::Integer,
    strains::Strains,
    behaviour::Behaviour;
    intervention::Intervention=DEFAULT_INTERVENTION,
    mode::Symbol=:SIR,
    arrivals::Integer=0
)::ParametricData where T <: AbstractRNG
    trials = length(rngs)

    susceptible = ParametricArray(zeros(
        Int, TIME_HORIZON + 1, trials, length(strains.initial), length(strains.strength),
        length(strains.radius), length(strains.infection_mean),
        length(strains.infection_shape)
    ))
    infected = ParametricArray(zeros(
        Int, TIME_HORIZON + 1, trials, length(strains.initial), length(strains.strength),
        length(strains.radius), length(strains.infection_mean),
        length(strains.infection_shape)
    ))
    recovered = ParametricArray(zeros(
        Int, TIME_HORIZON + 1, trials, length(strains.initial), length(strains.strength),
        length(strains.radius), length(strains.infection_mean),
        length(strains.infection_shape)
    ))

    for (i, initial) in enumerate(strains.initial),
            (j, strength) in enumerate(strains.strength),
            (k, radius) in enumerate(strains.radius),
            (h, mean) in enumerate(strains.infection_mean),
            (l, shape) in enumerate(strains.infection_shape)
        strain = Strain(initial, strength, radius, mean, shape)
        data = simulate(
            rngs, population, strain, behaviour, intervention=intervention, mode=mode,
            arrivals=arrivals
        )

        susceptible[
            initial=i, strength=j, radius=k, infection_mean=h, infection_shape=l
        ] = data.susceptible
        infected[
            initial=i, strength=j, radius=k, infection_mean=h, infection_shape=l
        ] = data.infected
        recovered[
            initial=i, strength=j, radius=k, infection_mean=h, infection_shape=l
        ] = data.recovered
    end

    return (
        population=population + arrivals,
        strains=strains,
        susceptible=susceptible,
        infected=infected,
        recovered=recovered
    )
end

function simulate(
    population::Integer,
    strain::Union{Strain, Strains},
    behaviour::Behaviour;
    kwargs...
)
    return simulate([GLOBAL_RNG], population, strain, behaviour; kwargs...)
end

function simulate(
    rngs::Vector{T},
    population::Integer,
    strain::Union{Strain, Strains};
    kwargs...
) where T <: AbstractRNG
    return simulate(rngs, population, strain, GLOBAL_BEHAVIOUR; kwargs...)
end

function simulate(population::Integer, strain::Union{Strain, Strains}; kwargs...)
    return simulate([GLOBAL_RNG], population, strain, GLOBAL_BEHAVIOUR; kwargs...)
end

function simulate(
    rng::AbstractRNG,
    population::Integer,
    strain::Union{Strain, Strains},
    behaviour::Behaviour;
    kwargs...
)
    return simulate([rng], population, strain, behaviour; kwargs...)
end

function simulate(
    rng::AbstractRNG,
    population::Integer,
    strain::Union{Strain, Strains};
    kwargs...
)
    return simulate([rng], population, strain, GLOBAL_BEHAVIOUR; kwargs...)
end
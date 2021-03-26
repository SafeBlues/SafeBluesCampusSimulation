using Random: AbstractRNG, Random.GLOBAL_RNG

using Distributions: Binomial, Discrete, Gamma, Sampleable, Univariate
using FileIO: load
using Images: Gray
using YAML: load_file

const HOURS_IN_DAY = 24
const DAYS_IN_WEEK = 7
hour(time::Integer)::Int = time % HOURS_IN_DAY + 1
day(time::Integer)::Int = (time ÷ HOURS_IN_DAY) % DAYS_IN_WEEK + 1
is_weekend(time::Integer)::Bool = day(time) == 6 || day(time) == 7

"""
    State(time, susceptible, infected, recovered, recovery_times)

Represents the epidemic state of a virus strain within a population.

**Fields**
- `time::Int`: The current time (measured in hours from 12:00am Monday).
- `susceptible::Int`: The number of susceptible individuals.
- `infected::Int`: The number of infected individuals.
- `recovered::Int`: The number of recovered individuals.
- `recovery_times::Vector{Float64}`: A collection (sorted in ascending order) of recovery
    times (measured in hours from 12:00am Monday).


    State(susceptible, infected, recovered)

Creates a `State` with the given number of susceptible, infected, and recovered individuals
and the default starting time `t = 0` and recovery times `recovery_times = []`.
"""
mutable struct State
    time::Int
    susceptible::Int
    infected::Int
    recovered::Int
    recovery_times::Vector{Float64}
end

function State(susceptible::Int, infected::Int, recovered::Int)
    return State(0, susceptible, infected, recovered, [])
end

"""
    Strain(strength, radius, recovery_shape, recovery_scale)

Stores the parameters describing the dynamics of a virus strain.

**Fields*
- `strength (Float64)`: The infection strength of the virus strain.
- `radius (Float64)`: The maximum infection radius of the virus strain (in metres).
- `recovery_shape (Float64)`: The shape parameter used by the gamma-distributed infection
    durations.
- `recovery_scale (Float64)`: The scale parameter used by the gamma-distributed infection
    durations.
"""
struct Strain
    strength::Float64
    radius::Float64
    recovery_shape::Float64
    recovery_scale::Float64
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
    index = Int(floor(u)) + 1

    if u - index + 1 >= sampler.probabilities[index]
        index = sampler.aliases[index]
    end

    x = sampler.scale * (((index - 1) % sampler.rows) + rand(rng))
    y = sampler.scale * (((index - 1) ÷ sampler.rows) + rand(rng))

    return (x=x, y=y)
end

"""
    Environment(campus_sampler, weekday_attendance, weekend_attendance, compliance)

Stores the population dynamics of a campus.

**Fields**
- `campus_sampler::CampusSampler`: Generates the positions of individuals on campus.
- `weekday_attendance::Vector{Float64}`: The attendance probabilities for each hour of a
    weekday.
- `weekend_attendance::Vector{Float64}`: The attendance probabilities for each hour of a
    weekend.
- `compliance::Float64`: The probability that an individual is running the Safe Blues app.
"""
struct Environment
    campus_sampler::CampusSampler
    weekday_attendance::Vector{Float64}
    weekend_attendance::Vector{Float64}
    compliance::Float64
end

"""
    load_environment(environment_file, heatmap_file)

Loads an `Environment` from an environment configuration and a campus heatmap.

**Arguments**
- `environment_file::String`: The path of an environment configuration file (`.yaml`).
- `heatmap_file::String`: The path of a campus heatmap file (`.png`).
"""
function load_environment(environment_file::String, heatmap_file::String)
    contents = load_file(environment_file)

    scale::Float64 = contents["scale"]
    weekday_attendance::Vector{Float64} = contents["weekday_attendance"]
    weekend_attendance::Vector{Float64} = contents["weekend_attendance"]
    compliance::Float64 = contents["compliance"]

    weights = (x -> Float64(Gray(x))).(load(heatmap_file))
    campus_sampler = CampusSampler(scale, weights)

    return Environment(campus_sampler, weekday_attendance, weekend_attendance, compliance)
end

const GLOBAL_ENVIRONMENT = cd(@__DIR__) do
    return load_environment("../environment.yaml", "../heatmap.png")
end

"""
    infect!(rng, state, strain, environment)

Handles the spread of a virus strain within the population.

Returns the number of new infections and updates `state` to reflect the new epidemic state
of the strain.

**Arguments**
- `rng::AbstractRNG`: A random number generator.
- `state::State`: The initial state of the strain.
- `strain::Strain`: Stores the parameters describing the virus strain.
- `environment::Environment`: Stores the population dynamics of the campus.
"""
function infect!(rng::AbstractRNG, state::State, strain::Strain, environment::Environment)
    # Get the number of active (attending & complying) participants.
    activity = (
        is_weekend(state.time) ? environment.weekend_attendance
        : environment.weekday_attendance
    )[hour(state.time)] * environment.compliance
    active_susceptible = rand(rng, Binomial(state.susceptible, activity))
    active_infected = rand(rng, Binomial(state.infected, activity))
    if active_susceptible == 0 || active_infected == 0
        return 0
    end

    width = environment.campus_sampler.columns * environment.campus_sampler.scale
    height = environment.campus_sampler.rows * environment.campus_sampler.scale
    blocks = zeros(Int, Int(height ÷ strain.radius) + 1, Int(width ÷ strain.radius) + 1)
    links = zeros(Int, active_infected)

    # Generate the locations of infected individuals.
    infected_points = [rand(rng, environment.campus_sampler) for _ in 1:active_infected]
    for (k, point) in enumerate(infected_points)
        i, j = Int(point.x ÷ strain.radius + 1), Int(point.y ÷ strain.radius + 1)
        blocks[i, j], links[k] = k, blocks[i, j]
    end

    bound = strain.radius^2

    offsets = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1))
    neighbours(i::Integer, j::Integer) = (map(+, (i, j), offset) for offset in offsets)

    # Attempt infections between nearby individuals.
    infections = 0
    for _ in 1:active_susceptible
        point = rand(rng, environment.campus_sampler)
        i, j = Int(point.x ÷ strain.radius + 1), Int(point.y ÷ strain.radius + 1)

        p = 0
        for (i′, j′) in neighbours(i, j)
            k = blocks[i′, j′]
            while k != 0
                point′ = infected_points[k]
                k = links[k]

                distance = (point.x - point′.x)^2 + (point.y - point′.y)^2
                if distance >= bound
                    continue
                end

                q = 1 - exp(-strain.strength * (1 - √distance / strain.radius))
                p = p + q - p * q
            end
        end

        # Determine whether an infection occurs.
        if rand(rng) < p
            infections += 1

            # Add the newly infected individual's recovery time to the recovery queue.
            t = state.time + rand(rng, Gamma(strain.recovery_shape, strain.recovery_scale))
            k = 1
            while k <= length(state.recovery_times)
                if t < state.recovery_times[k]
                    break
                end
                k += 1
            end
            insert!(state.recovery_times, k, t)
        end
    end

    state.susceptible -= infections
    state.infected += infections
    return infections
end

"""
    recover!(state)

Handles the recoveries within the population.

Returns the number of new recoveries and updates `state` to reflect the new epidemic state
of the strain.

**Arguments**
- `state::State`: The initial state of the strain.
"""
function recover!(state::State)
    recoveries = 0
    for t in state.recovery_times
        if t > state.time
            break
        end

        popfirst!(state.recovery_times)
        recoveries += 1
    end

    state.infected -= recoveries
    state.recovered += recoveries
    return recoveries
end

"""
    advance!(rng, state, strain, environment)

Advances the epidemic state of a virus strain forward by a single time increment.

Returns a `NamedTuple` storing the number of susceptible (`.susceptible`), infected
(`.infected`), and recovered (`.recovered`) individuals. Updates `state` to reflect the new
epidemic state of the virus.

**Arguments**
- `rng::AbstractRNG`: A random number generator.
- `state::State`: The initial state of the strain.
- `strain::Strain`: Stores the parameters describing the virus strain.
- `environment::Environment`: Stores the population dynamics of the campus.
"""
function advance!(rng::AbstractRNG, state::State, strain::Strain, environment::Environment)
    state.time += 1
    infect!(rng, state, strain, environment)
    recover!(state)

    return (
        susceptible=state.susceptible, infected=state.infected, recovered=state.recovered
    )
end

SimulationData = @NamedTuple begin
    population::Int
    susceptible::Matrix{Int}
    infected::Matrix{Int}
    recovered::Matrix{Int}
end

"""
    simulate(rngs, time_steps, state, strain, environment)

Simulates the spread of a virus strain within a population and returns `SimulationData`.

**Arguments**
- `rngs::Vector{AbstractRNG}`: A random number generator for each simulation trial.
- `time_steps::Integer`: The number of time steps within the simulation.
- `state::State`: The initial state of the strain.
- `strain::Strain`: Stores the parameters describing the virus strain.
- `environment::Environment`: Stores the population dynamics of the campus.


    simulate(time_steps, state, strain, environment)
    simulate(rngs, time_steps, state, strain)
    simulate(time_steps, state, strain)

Alternatively, when `rngs` is ommited the random number generators default to
`[Random.GLOBAL_RNG]` and when `environment` is ommited the environment defaults to
`GLOBAL_ENVIRONMENT`.


    simulate(rng, time_steps, state, strain, environment)
    simulate(rng, time_steps, state, strain)

If only a single random number generator `rng::AbstractRNG` is passed to `simulate`, then
the simulation will be run exactly once.
"""
function simulate(
    rngs::Vector{T},
    time_steps::Integer,
    state::State,
    strain::Strain,
    environment::Environment
)::SimulationData where T <: AbstractRNG
    trials = length(rngs)
    susceptible = zeros(Int, time_steps + 1, trials)
    infected = zeros(Int, time_steps + 1, trials)
    recovered = zeros(Int, time_steps + 1, trials)

    susceptible[1, :] = fill(state.susceptible, trials)
    infected[1, :] = fill(state.infected, trials)
    recovered[1, :] = fill(state.recovered, trials)

    for i in 1:trials
        trial_state = State(
            state.time, state.susceptible, state.infected, state.recovered,
            copy(state.recovery_times)
        )

        for j in 2:(time_steps + 1)
            data = advance!(rngs[i], trial_state, strain, environment)

            susceptible[j, i] = data.susceptible
            infected[j, i] = data.infected
            recovered[j, i] = data.recovered
        end
    end

    return (
        population=(state.susceptible + state.infected + state.recovered),
        susceptible=susceptible,
        infected=infected,
        recovered=recovered
    )
end

function simulate(
    time_steps::Integer,
    state::State,
    strain::Strain,
    environment::Environment
)
    return simulate([GLOBAL_RNG], time_steps, state, strain, environment)
end

function simulate(
    rngs::Vector{T},
    time_steps::Integer,
    state::State,
    strain::Strain
) where T <: AbstractRNG
    return simulate(rngs, time_steps, state, strain, GLOBAL_ENVIRONMENT)
end

function simulate(time_steps::Integer, state::State, strain::Strain)
    return simulate([GLOBAL_RNG], time_steps, state, strain, GLOBAL_ENVIRONMENT)
end

function simulate(
    rng::AbstractRNG,
    time_steps::Integer,
    state::State,
    strain::Strain,
    environment::Environment
)
    return simulate([rng], time_steps, state, strain, environment)
end

function simulate(rng::AbstractRNG, time_steps::Integer, state::State, strain::Strain)
    return simulate([rng], time_steps, state, strain, GLOBAL_ENVIRONMENT)
end
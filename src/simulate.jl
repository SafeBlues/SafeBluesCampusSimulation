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

struct Strain
    strength::Float64
    radius::Float64
    recovery_shape::Float64
    recovery_scale::Float64
end

struct CampusSampler <: Sampleable{Univariate, Discrete}
    width::Int
    height::Int
    scale::Float64
    probabilities::Matrix{Float64}
    aliases::Matrix{Int}
end

function CampusSampler(scale::Float64, weights::Matrix{Float64})
    width, height = size(weights)

    probabilities = weights / sum(weights) * width * height
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

    return CampusSampler(width, height, scale, probabilities, aliases)
end

function Base.rand(rng::AbstractRNG, s::CampusSampler)
    u = s.width * s.height * rand(rng)
    index = Int(floor(u)) + 1

    if u - index + 1 >= s.probabilities[index]
        index = s.aliases[index]
    end

    width, height = size(s.probabilities)

    x = s.scale * ((index % s.width) + rand(rng))
    y = s.scale * ((index ÷ s.height) + rand(rng))

    return (x=x, y=y)
end

struct Environment
    time_horizon::Int
    campus_sampler::CampusSampler
    weekday_attendance::Vector{Float64}
    weekend_attendance::Vector{Float64}
    compliance::Float64
end

function load_environment(environment_file::String, heatmap_file::String)
    contents = load_file(environment_file)

    time_horizon::Int = contents["time_horizon"]
    scale::Float64 = contents["scale"]
    weekday_attendance::Vector{Float64} = contents["weekday_attendance"]
    weekend_attendance::Vector{Float64} = contents["weekend_attendance"]
    compliance::Float64 = contents["compliance"]

    weights = (x -> Float64(Gray(x))).(load(heatmap_file))
    campus_sampler = CampusSampler(scale, weights)

    return Environment(
        time_horizon, campus_sampler, weekday_attendance, weekend_attendance, compliance
    )
end

const GLOBAL_ENVIRONMENT = cd(@__DIR__) do
    return load_environment("../environment.yaml", "../heatmap.png")
end

function infect!(rng::AbstractRNG, state::State, strain::Strain, environment::Environment)
    # Get the activity (attendance & compliance) probability.
    active_chance = (
        is_weekend(state.time) ? environment.weekend_attendance
        : environment.weekday_attendance
    )[hour(state.time)] * environment.compliance
    active_susceptible = rand(rng, Binomial(state.susceptible, active_chance))
    active_infected = rand(rng, Binomial(state.infected, active_chance))

    # Generate the location of every susceptible and infected individual.
    S = [rand(rng, environment.campus_sampler) for _ in 1:active_susceptible]
    I = [rand(rng, environment.campus_sampler) for _ in 1:active_infected]

    bound = strain.radius

    # Attempt infections between nearby individuals.
    infections = 0
    for s in S
        # Compute the probability of infection.
        p = 0
        for i in I
            distance = (s.x - i.x)^2 + (s.y - i.y)^2
            if distance >= bound
                continue
            end

            q = 1 - exp(-strain.strength * (1 - √distance / strain.radius))
            p = p + q - p * q
        end

        # Determine whether an infection occurs.
        if rand(rng) < p
            infections += 1

            # Add the newly infected person's recovery time to the recovery queue.
            t = state.time + rand(Gamma(strain.recovery_shape, strain.recovery_scale))
            i = 1
            while i <= length(state.recovery_times)
                if t < state.recovery_times[i]
                    break
                end
                i += 1
            end
            insert!(state.recovery_times, i, t)
        end
    end

    state.susceptible -= infections
    state.infected += infections
    return infections
end

function recover!(rng::AbstractRNG, state::State, strain::Strain)
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

function advance!(rng::AbstractRNG, state::State, strain::Strain, environment::Environment)
    state.time += 1
    infect!(rng, state, strain, environment)
    recover!(rng, state, strain)

    return (
        susceptible=state.susceptible, infected=state.infected, recovered=state.recovered
    )
end

function simulate!(rng::AbstractRNG, state::State, strain::Strain, environment::Environment)
    susceptible = zeros(Int, environment.time_horizon)
    infected = zeros(Int, environment.time_horizon)
    recovered = zeros(Int, environment.time_horizon)

    for i in 1:environment.time_horizon
        data = advance!(rng, state, strain, environment)

        susceptible[i] = data.susceptible
        infected[i] = data.infected
        recovered[i] = data.recovered
    end

    return (susceptible=susceptible, infected=infected, recovered=recovered)
end

function simulate!(state::State, strain::Strain, environment::Environment)
    return simulate!(GLOBAL_RNG, state, strain, environment)
end

function simulate!(rng::AbstractRNG, state::State, strain::Strain)
    return simulate!(rng, state, strain, GLOBAL_ENVIRONMENT)
end

function simulate!(state::State, strain::Strain)
    return simulate!(GLOBAL_RNG, state, strain, GLOBAL_ENVIRONMENT)
end
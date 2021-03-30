using NamedDims: unname
using PlotlyJS: AbstractTrace, Layout, Plot, attr, scatter, surface

const RED = "#FF595E"
const GREEN = "#8AC926"
const BLUE = "#1982C4"

const OPACITY = 0.1

"""
    sir_plot(data)

Plots an epidemic trajectory from the provided simulation data.

**Arguments**
- `data::SimulationData`: Stores the epidemic data generated by `simulate`.

**Keyword Arguments
- `show_susceptible::Bool=true`: Indicates whether the number of susceptible individuals is
    displayed.
- `show_infected::Bool=true`: Indicates whether the number of infected individuals is
    displayed.
- `show_recovered::Bool=true`: Indicates whether the number of recovered individuals is
    displayed.
- `show_trials::Bool=true`: Indicates whether the trajectories from individual trials are
    displayed.
"""
function sir_plot(
    data::SimulationData;
    show_susceptible::Bool=true,
    show_infected::Bool=true,
    show_recovered::Bool=true,
    show_trials::Bool=true
)
    traces::Vector{AbstractTrace} = []

    # Add the susceptible traces.
    times = (0:(size(data.susceptible, :time) - 1)) / HOURS_IN_DAY
    trials = size(data.susceptible, :trial)
    if show_susceptible && show_trials
        append!(traces, (scatter(
            x=times, y=data.susceptible[trial=trial], hoverinfo="skip", line_color=BLUE,
            mode="lines", opacity=OPACITY, showlegend=false
        ) for trial in 1:trials))
    end
    if show_susceptible
        push!(traces, scatter(
            x=times, y=vec(sum(data.susceptible, dims=:trial)) / trials, line_color=BLUE,
            mode="lines", name="Susceptible"
        ))
    end

    # Add the infected traces.
    times = (0:(size(data.infected, :time) - 1)) / HOURS_IN_DAY
    trials = size(data.infected, :trial)
    if show_infected && show_trials
        append!(traces, (scatter(
            x=times, y=data.infected[trial=trial], hoverinfo="skip", line_color=RED,
            mode="lines", opacity=OPACITY, showlegend=false
        ) for trial in 1:trials))
    end
    if show_infected
        push!(traces, scatter(
            x=times, y=vec(sum(data.infected, dims=:trial)) / trials, line_color=RED,
            mode="lines", name="Infected"
        ))
    end

    # Add the recovered traces.
    times = (0:(size(data.recovered, :time) - 1)) / HOURS_IN_DAY
    trials = size(data.recovered, :trial)
    if show_recovered && show_trials
        append!(traces, (scatter(
            x=times, y=data.recovered[trial=trial], hoverinfo="skip", line_color=GREEN,
            mode="lines", opacity=OPACITY, showlegend=false
        ) for trial in 1:trials))
    end
    if show_recovered
        push!(traces, scatter(
            x=times, y=vec(sum(data.recovered, dims=:trial)) / trials, line_color=GREEN,
            mode="lines", name="Recovered"
        ))
    end

    layout = Layout(
        xaxis_title="Time (days)", yaxis_range=(0, data.population),
        yaxis_title="Population"
    )

    return Plot(traces, layout)
end

"""
    cumulative_plot(data)

Plots the cumulative infections from the provided simulation data.

**Arguments**
- `data::SimulationData`: Stores the epidemic data generated by `simulate`.

**Keyword Arguments
- `show_trials::Bool=true`: Indicates whether the trajectories from individual trials are
    displayed.
"""
function cumulative_plot(data::SimulationData; show_trials::Bool=true)
    cumulative = data.population .- data.susceptible
    times = (0:(size(cumulative, :time) - 1)) / HOURS_IN_DAY
    trials = size(cumulative, :trial)

    traces::Vector{AbstractTrace} = []

    if show_trials
        append!(traces, (scatter(
            x=times, y=cumulative[trial=trial], hoverinfo="skip", line_color=RED,
            mode="lines", opacity=OPACITY, showlegend=false
        ) for trial in 1:trials))
    end

    push!(traces, scatter(
        x=times, y=vec(sum(cumulative, dims=:trial)) / trials, line_color=RED, mode="lines",
        name="", showlegend=false
    ))

    layout = Layout(
        xaxis_title="Time (days)", yaxis_range=(0, data.population),
        yaxis_title="Cumulative Infections"
    )

    return Plot(traces, layout)
end

"""
    reproduction_plot(data)

Plots the effective reproduction number from the provided simulation data.

**Arguments**
- `data::SimulationData`: Stores the epidemic data generated by `simulate`.

**Keyword Arguments
- `show_trials::Bool=true`: Indicates whether the trajectories from individual trials are
    displayed.
"""
function reproduction_plot(data::SimulationData, weight::Float64; show_trials::Bool=true)
    # Calculate the exponential moving average of the effective reproduction number.
    reproduction = copy(data.reproduction)
    for time in 1:(size(data.reproduction, :time) - 1)
        reproduction[time=(time + 1)] = (
            weight * reproduction[time=time] + (1 - weight) * reproduction[time=(time + 1)]
        )
    end

    times = (0:(size(data.reproduction, :time) - 1)) / HOURS_IN_DAY
    trials = size(data.reproduction, :trial)

    traces::Vector{AbstractTrace} = []

    # Add the traces for the effective reproduction number of individual trials.
    if show_trials
        append!(traces, (scatter(
            x=times, y=reproduction[trial=trial], hoverinfo="skip", line_color=BLUE,
            mode="lines", opacity=OPACITY, showlegend=false
        ) for trial in 1:trials))
    end

    # Add the average trace.
    push!(traces, scatter(
        x=times, y=vec(sum(reproduction, dims=:trial)) / trials, line_color=BLUE,
        mode="lines", name="", showlegend=false
    ))

    layout = Layout(
        xaxis_title="Time (days)", yaxis_title="Effective Reproduction Number"
    )

    return Plot(traces, layout)
end

const PARAMETER_LABELS = Dict(
    :strength => "Infection Strength",
    :radius => "Infection Radius (metres)",
    :shape => "Infection Duration Shape",
    :scale => "Infection Duration Scale"
)

"""
    parametric_plot(data, dim1, dim2)

Plots the total number of infections as a function of virus parameters.

**Arguments**
- `data::ParametricData`: Stores the epidemic data generated by `simulate`.
- `dim1::Symbol`: The parameter to show on the x-axis.
- `dim2::Symbol`: The parameter to show on the y-axis.

The dimensions `dim1` and `dim2` must be one of `:strength`, `:radius`, `:shape`, or
`:scale`. The keyword arguments `strength`, `radius`, `shape`, and `scale` control the
parameters that do not appear on the x-axis or y-axis.
"""
function parametric_plot(
    data::ParametricData,
    dim1::Symbol,
    dim2::Symbol;
    strength::Integer=1,
    radius::Integer=1,
    shape::Integer=1,
    scale::Integer=1
)
    # Get the indices of the relevant data.
    strength = :strength == dim1 || :strength == dim2 ? (:) : strength
    radius = :radius == dim1 || :radius == dim2 ? (:) : radius
    shape = :shape == dim1 || :shape == dim2 ? (:) : shape
    scale = :scale == dim1 || :scale == dim2 ? (:) : scale

    population = data.state.susceptible + data.state.infected + data.state.recovered
    total = population .- data.susceptible[
        time=end, strength=strength, radius=radius, shape=shape, scale=scale
    ]

    trials = size(total, :trial)
    parameters1 = getfield(data.strains, dim1)
    parameters2 = getfield(data.strains, dim2)

    # Add the average trace.
    traces = [surface(
        x=parameters1, y=parameters2, z=(sum(total, dims=:trial) / trials)[trial=1],
        colorscale=[[0.0, BLUE], [1.0, GREEN]]
    )]

    layout = Layout(scene=attr(
        xaxis_title=PARAMETER_LABELS[dim1], yaxis_title=PARAMETER_LABELS[dim2],
        zazis_range=(0, population), zaxis_title="Infected (Cumulative)"
    ))

    return Plot(traces, layout)
end
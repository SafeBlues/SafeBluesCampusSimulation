using Random: MersenneTwister

import Dash

using Dash: Input, Output, callback!, dash
using DashHtmlComponents
using DashCoreComponents
using YAML: load_file

config = cd(@__DIR__) do
    load_file("assets/config.yaml")
end

# ---------------------------------------------------------------------------------------- #
# Layout                                                                                   #
# ---------------------------------------------------------------------------------------- #
function make_marks(min::Real, max::Real)
    return Dict((i % 1.0 == 0 ? Int(i) : i) => "$i"
                for i in range(min, stop=max, length=5))
end

function app_header()
    return html_header(id="header", className="card card-header") do
        html_h1("Safe Blues Campus Simulation")
    end
end

function main_controls()
    return html_div(id="main-controls-card", className="card card-grid") do
        # Model Choice
        html_p(className="label", config["model_label"]),
        dcc_radioitems(
            id="model-input", className="radio-items", value=config["model_value"],
            options=config["model_options"]
        ),
        
        # Seed
        html_p(className="label", config["seed_label"]),
        dcc_input(
            id="seed-input", className="input-box", type="number",
            min=config["seed_min"], step=config["seed_step"], value=config["seed_value"]    
        ),

        # Trials
        html_p(className="label", config["trials_label"]),
        dcc_input(
            id="trials-input", className="input-box", type="number",
            min=config["trials_min"], step=config["trials_step"],
            value=config["trials_value"]
        ),

        # Population
        html_p(className="label", config["population_label"]),
        dcc_input(
            id="population-input", className="input-box", type="number",
            min=config["population_min"], step=config["population_step"],
            value=config["population_value"]
        ),

        # Initial Infection Chance
        html_p(className="label", config["infection_initial_label"]),
        dcc_input(
            id="infection-initial-input", className="input-box", type="number",
            min=config["infection_initial_min"], max=config["infection_initial_max"],
            value=config["infection_initial_value"]
        ),
        dcc_slider(
            id="infection-initial-slider", className="slider",
            min=config["infection_initial_min"], max=config["infection_initial_max"],
            step=config["infection_initial_step"], value=config["infection_initial_value"],
            marks=make_marks(
                config["infection_initial_min"], config["infection_initial_max"]
            )
        ),

        # Infection Strength
        html_p(className="label", config["infection_strength_label"]),
        dcc_input(
            id="infection-strength-input", className="input-box", type="number",
            min=config["infection_strength_min"], value=config["infection_strength_value"]
        ),
        dcc_slider(
            id="infection-strength-slider", className="slider",
            min=config["infection_strength_min"], max=config["infection_strength_max"],
            step=config["infection_strength_step"],
            value=config["infection_strength_value"], marks=make_marks(
                config["infection_strength_min"], config["infection_strength_max"]
            )
        ),

        # Infection Radius
        html_p(className="label", config["infection_radius_label"]),
        dcc_input(
            id="infection-radius-input", className="input-box", type="number",
            min=config["infection_radius_min"], value=config["infection_radius_value"]
        ),
        dcc_slider(
            id="infection-radius-slider", className="slider",
            min=config["infection_radius_min"], max=config["infection_radius_max"],
            step=config["infection_radius_step"],
            value=config["infection_radius_value"],
            marks=make_marks(config["infection_radius_min"], config["infection_radius_max"])
        ),

        # Infection Duration Mean
        html_p(className="label", config["infection_mean_label"]),
        dcc_input(
            id="infection-mean-input", className="input-box", type="number",
            min=config["infection_mean_min"], value=config["infection_mean_value"]
        ),
        dcc_slider(
            id="infection-mean-slider", className="slider",
            min=config["infection_mean_min"], max=config["infection_mean_max"],
            step=config["infection_mean_step"],
            value=config["infection_mean_value"],
            marks=make_marks(config["infection_mean_min"], config["infection_mean_max"])
        ),

        # Infection Duration Shape
        html_p(className="label", config["infection_shape_label"]),
        dcc_input(
            id="infection-shape-input", className="input-box", type="number",
            min=config["infection_shape_min"], value=config["infection_shape_value"]
        ),
        dcc_slider(
            id="infection-shape-slider", className="slider",
            min=config["infection_shape_min"], max=config["infection_shape_max"],
            step=config["infection_shape_step"],
            value=config["infection_shape_value"],
            marks=make_marks(config["infection_shape_min"], config["infection_shape_max"])
        ),

        # Intervention Start
        html_p(className="label", config["intervention_start_label"]),
        dcc_input(
            id="intervention-start-input", className="input-box", type="number",
            min=config["intervention_start_min"], max=config["intervention_start_max"],
            step=config["intervention_start_step"],
            value=config["intervention_start_value"]
        ),

        # Intervention Stop
        html_p(className="label", config["intervention_stop_label"]),
        dcc_input(
            id="intervention-stop-input", className="input-box", type="number",
            min=config["intervention_stop_min"], max=config["intervention_stop_max"],
            step=config["intervention_stop_step"],
            value=config["intervention_stop_value"]
        ),

        # Intervention Strength Scale
        html_p(className="label", config["intervention_strength_label"]),
        dcc_input(
            id="intervention-strength-input", className="input-box", type="number",
            min=config["intervention_strength_min"],
            value=config["intervention_strength_value"]
        ),
        dcc_slider(
            id="intervention-strength-slider", className="slider",
            min=config["intervention_strength_min"],
            max=config["intervention_strength_max"],
             step=config["intervention_strength_step"],
            value=config["intervention_strength_value"],
            marks=make_marks(
                config["intervention_strength_min"], config["intervention_strength_max"]
            )
        )
    end
end

function parametric_controls()
    return html_div(id="parametric-controls-card", className="card card-grid") do
        # Parameter Choice
        html_p(className="label", config["parameter_label"]),
        dcc_dropdown(
            id="parameter-input", value=config["parameter_value"],
            options=config["parameter_options"], multi=true
        ),

        # Resolution
        html_p(className="label", config["resolution_label"]),
        dcc_input(
            id="parametric-resolution-input", className="input-box", type="number",
            min=config["resolution_min"], step=config["resolution_step"],
            value=config["resolution_value"]
        ),

        # Trials
        html_p(className="label", config["trials_label"]),
        dcc_input(
            id="parametric-trials-input", className="input-box", type="number",
            min=config["trials_min"], step=config["trials_step"],
            value=config["trials_value"]
        ),

        # Infection Strength
        html_p(className="label", config["infection_strength_label"]),
        dcc_input(
            id="infection-strength-start", className="input-box-small", type="number",
            min=config["infection_strength_min"],
            value=config["infection_strength_start_value"]
        ),
        dcc_input(
            id="infection-strength-stop", className="input-box-small", type="number",
            min=config["infection_strength_min"],
            value=config["infection_strength_stop_value"]
        ),
        dcc_rangeslider(
            id="infection-strength-range", className="slider",
            min=config["infection_strength_min"], max=config["infection_strength_max"],
            step=config["infection_strength_step"], value=(
                config["infection_strength_start_value"],
                config["infection_strength_stop_value"]
            ), marks=make_marks(
                config["infection_strength_min"], config["infection_strength_max"]
            )
        ),

        # Infection Radius
        html_p(className="label", config["infection_radius_label"]),
        dcc_input(
            id="infection-radius-start", className="input-box-small", type="number",
            min=config["infection_radius_min"],
            value=config["infection_radius_start_value"]
        ),
        dcc_input(
            id="infection-radius-stop", className="input-box-small", type="number",
            min=config["infection_radius_min"],
            value=config["infection_radius_stop_value"]
        ),
        dcc_rangeslider(
            id="infection-radius-range", className="slider",
            min=config["infection_radius_min"], max=config["infection_radius_max"],
            step=config["infection_radius_step"], value=(
                config["infection_radius_start_value"],
                config["infection_radius_stop_value"]
            ), marks=make_marks(
                config["infection_radius_min"], config["infection_radius_max"]
            )
        ),

        # Infection Duration Mean
        html_p(className="label", config["infection_mean_label"]),
        dcc_input(
            id="infection-mean-start", className="input-box-small", type="number",
            min=config["infection_mean_min"],
            value=config["infection_mean_start_value"]
        ),
        dcc_input(
            id="infection-mean-stop", className="input-box-small", type="number",
            min=config["infection_mean_min"],
            value=config["infection_mean_stop_value"]
        ),
        dcc_rangeslider(
            id="infection-mean-range", className="slider",
            min=config["infection_mean_min"], max=config["infection_mean_max"],
            step=config["infection_mean_step"], value=(
                config["infection_mean_start_value"],
                config["infection_mean_stop_value"]
            ), marks=make_marks(
                config["infection_mean_min"], config["infection_mean_max"]
            )
        ),

        # Infection Duration Shape
        html_p(className="label", config["infection_shape_label"]),
        dcc_input(
            id="infection-shape-start", className="input-box-small", type="number",
            min=config["infection_shape_min"],
            value=config["infection_shape_start_value"]
        ),
        dcc_input(
            id="infection-shape-stop", className="input-box-small", type="number",
            min=config["infection_shape_min"],
            value=config["infection_shape_stop_value"]
        ),
        dcc_rangeslider(
            id="infection-shape-range", className="slider",
            min=config["infection_shape_min"], max=config["infection_shape_max"],
            step=config["infection_shape_step"], value=(
                config["infection_shape_start_value"],
                config["infection_shape_stop_value"]
            ), marks=make_marks(
                config["infection_shape_min"], config["infection_shape_max"]
            )
        )
    end
end

function draw_sir_plot()
    return html_div(id="sir-plot-card", className="card") do
        dcc_graph(id="sir-plot")
    end
end

function draw_cumulative_plot()
    return html_div(id="cumulative-plot-card", className="card") do
        dcc_graph(id="cumulative-plot")
    end
end

function draw_parametric_plot()
    return html_div(id="parametric-plot-card", className="card") do
        dcc_graph(id="parametric-plot")
    end
end

function app_layout()
    return html_div(className="app-grid") do
        app_header(),
        main_controls(),
        draw_sir_plot(),
        draw_cumulative_plot(),
        draw_parametric_plot(),
        parametric_controls()
    end
end

app = cd(@__DIR__) do
    return dash()
end

app.layout = app_layout()

# ---------------------------------------------------------------------------------------- #
# Callbacks                                                                                #
# ---------------------------------------------------------------------------------------- #
function sync_slider(input_id::String, slider_id::String)
    # Send the current value from the slider to the input box.
    callback!(
        app, Output(input_id, "value"), Input(slider_id, "value")
    ) do value
        return value
    end

    # Send the current value from the input box to the slider.
    callback!(
        app, Output(slider_id, "value"), Input(input_id, "value")
    ) do value
        return value
    end
end

function sync_rangeslider(input1_id::String, input2_id::String, slider_id::String)
    # Send the current value from the slider to the input boxes.
    callback!(
        app,
        Output(input1_id, "value"),
        Output(input2_id, "value"),
        Input(slider_id, "value")
    ) do value
        return value
    end

    # Send the current value from the input boxes to the slider.
    callback!(
        app,
        Output(slider_id, "value"),
        Input(input1_id, "value"),
        Input(input2_id, "value")
    ) do value1, value2
        return (value1, value2)
    end
end

# Connect the parameter inputs to the parameter sliders.
sync_slider("infection-initial-input", "infection-initial-slider")
sync_slider("infection-strength-input", "infection-strength-slider")
sync_slider("infection-radius-input", "infection-radius-slider")
sync_slider("infection-mean-input", "infection-mean-slider")
sync_slider("infection-shape-input", "infection-shape-slider")
sync_slider("intervention-strength-input", "intervention-strength-slider")

sync_rangeslider(
    "infection-strength-start", "infection-strength-stop", "infection-strength-range"
)
sync_rangeslider(
    "infection-radius-start", "infection-radius-stop", "infection-radius-range"
)
sync_rangeslider(
    "infection-mean-start", "infection-mean-stop", "infection-mean-range"
)
sync_rangeslider(
    "infection-shape-start", "infection-shape-stop", "infection-shape-range"
)

# Disable the infection duration sliders when "SI" model is used.
callback!(
    app,
    Output("infection-mean-input", "disabled"),
    Output("infection-mean-slider", "disabled"),
    Output("infection-shape-input", "disabled"),
    Output("infection-shape-slider", "disabled"),
    Input("model-input", "value")
) do model
    disabled = model != "SIR"
    return (disabled, disabled, disabled, disabled)
end

# Connect the parameter inputs to the plots.
callback!(
    app,
    Output("sir-plot", "figure"),
    Output("cumulative-plot", "figure"),
    Input("model-input", "value"),
    Input("seed-input", "value"),
    Input("trials-input", "value"),
    Input("population-input", "value"),
    Input("infection-initial-input", "value"),
    Input("infection-strength-input", "value"),
    Input("infection-radius-input", "value"),
    Input("infection-mean-input", "value"),
    Input("infection-shape-input", "value"),
    Input("intervention-start-input", "value"),
    Input("intervention-stop-input", "value"),
    Input("intervention-strength-input", "value")
) do model, seed, trials, population, infection_initial, infection_strength,
        infection_radius, infection_mean, infection_shape, intervention_start,
        intervention_stop,
        intervention_strength
    if model == "SI"
        infection_shape = 1
        infection_mean = Inf
    end

    rngs = [MersenneTwister(seed + (i - 1) * trials) for i in 1:trials]
    strain = Strain(
        infection_initial, infection_strength, infection_radius, infection_mean,
        infection_shape
    )
    intervention = Intervention(intervention_start, intervention_stop, intervention_strength)

    data = simulate(rngs, population, strain; intervention=intervention)

    return sir_plot(data; show_recovered=(model == "SIR")), cumulative_plot(data)
end

# Enable the appropriate sliders when a parameter is selected.
callback!(
    app,
    Output("infection-strength-start", "disabled"),
    Output("infection-strength-stop", "disabled"),
    Output("infection-strength-range", "disabled"),
    Input("parameter-input", "value")
) do value
    if value == "" || !("strength" in value)
        return (true, true, true)
    end

    return (false, false, false)
end

callback!(
    app,
    Output("infection-radius-start", "disabled"),
    Output("infection-radius-stop", "disabled"),
    Output("infection-radius-range", "disabled"),
    Input("parameter-input", "value")
) do value
    if value == "" || !("radius" in value)
        return (true, true, true)
    end

    return (false, false, false)
end

callback!(
    app,
    Output("infection-mean-start", "disabled"),
    Output("infection-mean-stop", "disabled"),
    Output("infection-mean-range", "disabled"),
    Input("parameter-input", "value")
) do value
    if value == "" || !("duration-mean" in value)
        return (true, true, true)
    end

    return (false, false, false)
end

callback!(
    app,
    Output("infection-shape-start", "disabled"),
    Output("infection-shape-stop", "disabled"),
    Output("infection-shape-range", "disabled"),
    Input("parameter-input", "value")
) do value
    if value == "" || !("duration-shape" in value)
        return (true, true, true)
    end

    return (false, false, false)
end

# Connect the parameter inputs to the parametric plot.
callback!(
    app,
    Output("parametric-plot", "figure"),
    Input("model-input", "value"),
    Input("seed-input", "value"),
    Input("parametric-resolution-input", "value"),
    Input("parametric-trials-input", "value"),
    Input("population-input", "value"),
    Input("infection-initial-input", "value"),
    Input("infection-strength-input", "value"),
    Input("infection-strength-start", "value"),
    Input("infection-strength-stop", "value"),
    Input("infection-radius-input", "value"),
    Input("infection-radius-start", "value"),
    Input("infection-radius-stop", "value"),
    Input("infection-mean-input", "value"),
    Input("infection-mean-start", "value"),
    Input("infection-mean-stop", "value"),
    Input("infection-shape-input", "value"),
    Input("infection-shape-start", "value"),
    Input("infection-shape-stop", "value"),
    Input("intervention-start-input", "value"),
    Input("intervention-stop-input", "value"),
    Input("intervention-strength-input", "value"),
    Input("parameter-input", "value")
) do model, seed, resolution, trials, population, initial, infection_strength,
        infection_strength_start, infection_strength_stop, infection_radius,
        infection_radius_start, infection_radius_stop, infection_mean,
        infection_mean_start, infection_mean_stop, infection_shape, infection_shape_start,
        infection_shape_stop, intervention_start, intervention_stop, intervention_strength,
        parameters
    if parameters == "" || (length(parameters) != 1 && length(parameters) != 2)
        return Plot()
    end

    rngs = [MersenneTwister(seed + (i - 1) * trials) for i in 1:trials]
    intervention = Intervention(intervention_start, intervention_stop, intervention_strength)

    if "strength" in parameters
        infection_strength = Array(range(
            infection_strength_start, stop=infection_strength_stop, length=resolution
        ))
    end

    if "radius" in parameters
        infection_radius = Array(range(
            infection_radius_start, stop=infection_radius_stop, length=resolution
        ))
    end

    if "duration-mean" in parameters
        infection_mean = Array(range(
            infection_mean_start, stop=infection_mean_stop, length=resolution
        ))
    end

    if "duration-shape" in parameters
        infection_shape = Array(range(
            infection_shape_start, stop=infection_shape_stop, length=resolution
        ))
    end

    strains = Strains(
        initial, infection_strength, infection_radius, infection_mean, infection_shape    
    )

    data = simulate(rngs, population, strains; intervention=intervention)

    symbols = Dict(
        "strength" => :strength,
        "radius" => :radius,
        "duration-mean" => :duration_mean,
        "duration-shape" => :duration_shape
    )
    dims = map(k -> symbols[k], parameters)

    return parametric_plot(data, dims...)
end
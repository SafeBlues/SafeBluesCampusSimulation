using Random: MersenneTwister

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

function controls()
    return html_div(id="controls-card", className="card card-grid") do
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

function app_layout()
    return html_div(className="app-grid") do
        app_header(),
        controls(),
        draw_sir_plot(),
        draw_cumulative_plot()
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

# Connect the parameter inputs to the parameter sliders.
sync_slider("infection-initial-input", "infection-initial-slider")
sync_slider("infection-strength-input", "infection-strength-slider")
sync_slider("infection-radius-input", "infection-radius-slider")
sync_slider("infection-mean-input", "infection-mean-slider")
sync_slider("infection-shape-input", "infection-shape-slider")
sync_slider("intervention-strength-input", "intervention-strength-slider")

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
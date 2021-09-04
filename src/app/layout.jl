using DashHtmlComponents
using DashCoreComponents
using YAML: load_file

config = cd(@__DIR__) do
    load_file("config.yaml")
end

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
        dcc_dropdown(
            id="model-input", value=config["model_value"], options=config["model_options"]
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

        # Arrivals
        html_p(className="label", config["arrivals_label"]),
        dcc_input(
            id="arrivals-input", className="input-box", type="number",
            min=config["arrivals_min"], step=config["arrivals_step"],
            value=config["arrivals_value"]
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

        # Incubation Duration Mean
        html_p(className="label", config["incubation_mean_label"]),
        dcc_input(
            id="incubation-mean-input", className="input-box", type="number",
            min=config["incubation_mean_min"], value=config["incubation_mean_value"]
        ),
        dcc_slider(
            id="incubation-mean-slider", className="slider",
            min=config["incubation_mean_min"], max=config["incubation_mean_max"],
            step=config["incubation_mean_step"],
            value=config["incubation_mean_value"],
            marks=make_marks(config["incubation_mean_min"], config["incubation_mean_max"])
        ),

        # Incubation Duration Shape
        html_p(className="label", config["incubation_shape_label"]),
        dcc_input(
            id="incubation-shape-input", className="input-box", type="number",
            min=config["incubation_shape_min"], value=config["incubation_shape_value"]
        ),
        dcc_slider(
            id="incubation-shape-slider", className="slider",
            min=config["incubation_shape_min"], max=config["incubation_shape_max"],
            step=config["incubation_shape_step"],
            value=config["incubation_shape_value"],
            marks=make_marks(config["incubation_shape_min"], config["incubation_shape_max"])
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

function draw_infection_probability_plot()
    return html_div(id="infection-probability-plot-card", className="card") do
        dcc_graph(id="infection-probability-plot")
    end
end

function draw_duration_plot()
    return html_div(id="duration-plot-card", className="card") do
        dcc_graph(id="duration-plot")
    end
end

function app_layout()
    return html_div(className="app-grid") do
        app_header(),
        main_controls(),
        draw_infection_probability_plot(),
        draw_duration_plot(),
        draw_sir_plot()
    end
end

app.layout = app_layout()
using DashHtmlComponents
using DashCoreComponents

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
        html_p(className="label", "Model"),
        dcc_dropdown(
            id="model-input",
            value="SEIR",
            options=[
                Dict("label" => "SEIR", "value" => "SEIR"),
                Dict("label" => "SIR", "value" => "SIR"),
                Dict("label" => "SEIS", "value" => "SEIS"),
                Dict("label" => "SIS", "value" => "SIS"),
                Dict("label" => "SEI", "value" => "SEI"),
                Dict("label" => "SI", "value" => "SI")
            ]
        ),

        # Seed
        html_p(className="label", "Seed"),
        dcc_input(
            id="seed-input",
            className="input-box",
            type="number",
            min=0,
            step=1,
            value=0  
        ),

        # Trials
        html_p(className="label", "Trials"),
        dcc_input(
            id="trials-input",
            className="input-box",
            type="number",
            min=1,
            max=50,
            step=1,
            value=1
        ),

        # Population
        html_p(className="label", "Population"),
        dcc_input(
            id="population-input",
            className="input-box",
            type="number",
            min=1,
            max=500,
            step=1,
            value=250
        ),

        # Initial Infection Chance
        html_p(className="label", "Initial Infection"),
        dcc_input(
            id="infection-initial-input",
            className="input-box",
            type="number",
            min=0,
            max=1,
            value=0.1
        ),
        dcc_slider(
            id="infection-initial-slider",
            className="slider",
            min=0.0,
            max=1.0,
            step=0.01,
            value=0.1,
            marks=make_marks(0.0, 1.0)
        ),

        # Infection Strength
        html_p(className="label", "Infection Strength"),
        dcc_input(
            id="infection-strength-input",
            className="input-box",
            type="number",
            min=0.0, value=0.05
        ),
        dcc_slider(
            id="infection-strength-slider",
            className="slider",
            min=0.0,
            max=0.1,
            step=0.01,
            value=0.05,
            marks=make_marks(0.0, 0.1)
        ),

        # Infection Radius
        html_p(className="label", "Infection Radius"),
        dcc_input(
            id="infection-radius-input",
            className="input-box",
            type="number",
            min=0.0,
            value=10.0
        ),
        dcc_slider(
            id="infection-radius-slider",
            className="slider",
            min=0.0,
            max=20.0,
            step=0.01,
            value=10.0,
            marks=make_marks(0.0, 20.0)
        ),

        # Incubation Duration Mean
        html_p(className="label", "Incubation Mean"),
        dcc_input(
            id="incubation-mean-input",
            className="input-box",
            type="number",
            min=0.0,
            value=48.0
        ),
        dcc_slider(
            id="incubation-mean-slider",
            className="slider",
            min=0.0,
            max=336.0,
            step=0.01,
            value=48.0,
            marks=make_marks(0.0, 336.0)
        ),

        # Incubation Duration Shape
        html_p(className="label", "Incubation Shape"),
        dcc_input(
            id="incubation-shape-input",
            className="input-box",
            type="number",
            min=0.0,
            value=10.0
        ),
        dcc_slider(
            id="incubation-shape-slider",
            className="slider",
            min=0.0,
            max=20.0,
            step=0.01,
            value=10.0,
            marks=make_marks(0.0, 20.0)
        ),

        # Infection Duration Mean
        html_p(className="label", "Infection Mean"),
        dcc_input(
            id="infection-mean-input",
            className="input-box",
            type="number",
            min=0.0,
            value=168.0
        ),
        dcc_slider(
            id="infection-mean-slider",
            className="slider",
            min=0.0,
            max=336.0,
            step=0.01,
            value=168.0,
            marks=make_marks(0.0, 336.0)
        ),

        # Infection Duration Shape
        html_p(className="label", "Infection Shape"),
        dcc_input(
            id="infection-shape-input",
            className="input-box",
            type="number",
            min=0.0,
            value=10.0
        ),
        dcc_slider(
            id="infection-shape-slider",
            className="slider",
            min=0.0,
            max=20.0,
            step=0.01,
            value=10.0,
            marks=make_marks(0.0, 20.0)
        ),

        # Intervention Start
        html_p(className="label", "Intervention Start"),
        dcc_input(
            id="intervention-start-input",
            className="input-box",
            type="number",
            min=0,
            max=540,
            step=1,
            value=0
        ),

        # Intervention Stop
        html_p(className="label", "Intervention Stop"),
        dcc_input(
            id="intervention-stop-input",
            className="input-box",
            type="number",
            min=0,
            max=540,
            step=1,
            value=0
        ),

        # Intervention Strength Scale
        html_p(className="label", "Intervention Strength"),
        dcc_input(
            id="intervention-strength-input",
            className="input-box",
            type="number",
            min=0.0,
            max=1.0,
            value=0.5
        ),
        dcc_slider(
            id="intervention-strength-slider", className="slider",
            min=0.0,
            max=1.0,
            step=0.01,
            value=0.5,
            marks=make_marks(0.0, 1.0)
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
using Random: MersenneTwister

using Dash: Input, Output, callback!
using DashHtmlComponents
using DashCoreComponents

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
sync_slider("incubation-mean-input", "incubation-mean-slider")
sync_slider("incubation-shape-input", "incubation-shape-slider")
sync_slider("infection-mean-input", "infection-mean-slider")
sync_slider("infection-shape-input", "infection-shape-slider")
sync_slider("intervention-strength-input", "intervention-strength-slider")

# Disable the incubation duration sliders when "SIR", "SIS" or "SI" model is used.
callback!(
    app,
    Output("incubation-mean-input", "disabled"),
    Output("incubation-mean-slider", "disabled"),
    Output("incubation-shape-input", "disabled"),
    Output("incubation-shape-slider", "disabled"),
    Input("model-input", "value")
) do model
    disabled = model == "SIR" || model == "SIS" || model == "SI"
    return (disabled, disabled, disabled, disabled)
end


# Disable the infection duration sliders when "SI" or "SEI" model is used.
callback!(
    app,
    Output("infection-mean-input", "disabled"),
    Output("infection-mean-slider", "disabled"),
    Output("infection-shape-input", "disabled"),
    Output("infection-shape-slider", "disabled"),
    Input("model-input", "value")
) do model
    disabled = model == "SI" || model == "SEI"
    return (disabled, disabled, disabled, disabled)
end

# Connect the parameter inputs to the plots.
callback!(
    app,
    Output("sir-plot", "figure"),
    Output("infection-probability-plot", "figure"),
    Output("duration-plot", "figure"),
    Input("model-input", "value"),
    Input("seed-input", "value"),
    Input("trials-input", "value"),
    Input("population-input", "value"),
    Input("arrivals-input", "value"),
    Input("infection-initial-input", "value"),
    Input("infection-strength-input", "value"),
    Input("infection-radius-input", "value"),
    Input("incubation-mean-input", "value"),
    Input("incubation-shape-input", "value"),
    Input("infection-mean-input", "value"),
    Input("infection-shape-input", "value"),
    Input("intervention-start-input", "value"),
    Input("intervention-stop-input", "value"),
    Input("intervention-strength-input", "value")
) do model, seed, trials, population, arrivals, infection_initial, infection_strength,
        infection_radius, incubation_mean, incubation_shape, infection_mean,
        infection_shape, intervention_start, intervention_stop, intervention_strength
    if model == "SI"
        infection_shape = 1
        infection_mean = Inf
        model = :SI
    elseif model == "SIS"
        model = :SIS
    elseif model == "SIR"
        model = :SIR
    elseif model == "SEI"
        infection_shape = 1
        infection_mean = Inf
        model = :SEI
    elseif model == "SEIS"
        model = :SEIS
    elseif model == "SEIR"
        model = :SEIR
    end

    rngs = [MersenneTwister(seed + (i - 1) * trials) for i in 1:trials]
    strain = Strain(
        infection_initial, infection_strength, infection_radius, incubation_mean,
        incubation_shape, infection_mean, infection_shape
    )
    intervention = Intervention(intervention_start, intervention_stop, intervention_strength)

    data = simulate(
        rngs, population, strain; intervention=intervention, mode=model, arrivals=arrivals
    )

    return (
        sir_plot(
            data;
            show_exposed=(model == :SEIR || model == :SEI || model == :SEIS),
            show_recovered=(model == :SIR || model == :SEIR)
        ),
        infection_probability_plot(strain),
        duration_plot(
            strain;
            show_incubation=(model == :SEIR || model == :SEI || model == :SEIS)
        )
    )
end
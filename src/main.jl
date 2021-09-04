using Dash: dash, run_server

include("core/simulate.jl")
include("core/plots.jl")

if isinteractive()
    app = cd(@__DIR__) do
        return dash()
    end
else
    app = dash()
end

# App Layout
include("app/layout.jl")

# App Callbacks
include("app/callbacks.jl")

run_server(app, "0.0.0.0")
using Dash: run_server

include("core/simulate.jl")
include("core/plots.jl")

include("app/app.jl")

run_server(app, "0.0.0.0")
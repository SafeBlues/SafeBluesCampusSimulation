using SafeBluesCampusSimulation
using Documenter

DocMeta.setdocmeta!(SafeBluesCampusSimulation, :DocTestSetup, :(using SafeBluesCampusSimulation); recursive=true)

makedocs(;
    modules=[SafeBluesCampusSimulation],
    authors="Thomas Graham",
    repo="https://github.com/tjgraham/SafeBluesCampusSimulation.jl/blob/{commit}{path}#{line}",
    sitename="SafeBluesCampusSimulation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

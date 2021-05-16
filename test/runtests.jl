using Test

cd(@__DIR__) do
    include("../src/core/simulate.jl")

    @testset "SafeBluesCampusSimulation" begin
        include("testsets/time.jl")
        include("testsets/campus_sampler.jl")
    end
end

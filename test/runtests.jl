using SafeBluesCampusSimulation
using Test

cd(@__DIR__) do
    @testset "SafeBluesCampusSimulation.jl" begin
        include("testsets/time.jl")
        include("testsets/campus_sampler.jl")
    end
end

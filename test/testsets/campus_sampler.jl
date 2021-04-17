# Tests for `CampusSampler` defined in `src/core/simulate.jl`.

weights1 = [
    1.0 2.0;
    3.0 4.0;
]
probabilities1 = weights1 / sum(weights1)
sampler1 = SafeBluesCampusSimulation.CampusSampler(1.0, weights1)

weights2 = [
    0.0 1.0 2.0 1.0 0.0;
    1.0 2.0 3.0 2.0 1.0;
    2.0 3.0 4.0 3.0 2.0;
    1.0 2.0 3.0 2.0 1.0;
    0.0 1.0 2.0 1.0 0.0;
]
probabilities2 = weights2 / sum(weights2)
sampler2 = SafeBluesCampusSimulation.CampusSampler(3.0, weights2)

weights3 = [
    10.0 20.0 10.0 30.0 10.0 20.0 10.0
    20.0 30.0 35.0 40.0 35.0 30.0 20.0
    10.0 20.0 10.0 30.0 10.0 20.0 10.0
]
probabilities3 = weights3 / sum(weights3)
sampler3 = SafeBluesCampusSimulation.CampusSampler(10.0, weights3)

function compute_probabilities(sampler::SafeBluesCampusSimulation.CampusSampler)
    indices = LinearIndices((sampler.rows, sampler.columns))

    return (i -> sampler.probabilities[i] + sum(Float64[
        1 - sampler.probabilities[j] for j in indices if sampler.aliases[j] == i
    ])).(indices) / (sampler.rows * sampler.columns)
end

function estimate_probabilities(
    sampler::SafeBluesCampusSimulation.CampusSampler;
    N::Integer=1000000
)
    count = zeros(Int, sampler.rows, sampler.columns)
    for _ in 1:N
        i, j = (x -> Int(floor(x / sampler.scale)) + 1).(Tuple(rand(sampler)))
        count[i, j] += 1
    end

    return count / N
end

@testset "CampusSampler" begin
    @test sampler1.rows == 2
    @test sampler1.columns == 2
    @test sampler1.scale == 1.0
    @test compute_probabilities(sampler1) ≈ probabilities1
    @test isapprox(estimate_probabilities(sampler1), probabilities1, rtol=0.01)

    @test sampler2.rows == 5
    @test sampler2.columns == 5
    @test sampler2.scale == 3.0
    @test compute_probabilities(sampler2) ≈ probabilities2
    @test isapprox(estimate_probabilities(sampler2), probabilities2, rtol=0.01)

    @test sampler3.rows == 3
    @test sampler3.columns == 7
    @test sampler3.scale == 10.0
    @test compute_probabilities(sampler3) ≈ probabilities3
    @test isapprox(estimate_probabilities(sampler3), probabilities3, rtol=0.01)
end
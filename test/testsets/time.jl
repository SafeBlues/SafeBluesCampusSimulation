# Tests for time-related functions defined in `src/core/simulate.jl`.

@testset "Time" begin
    @test SafeBluesCampusSimulation.hour(1) == 1    # 1 -> 12am -- 1am
    @test SafeBluesCampusSimulation.hour(25) == 1   # 25 -> 12am -- 1am
    @test SafeBluesCampusSimulation.hour(50) == 2   # 50 -> 1am -- 2am
    @test SafeBluesCampusSimulation.hour(100) == 4  # 100 -> 3am -- 4am
    @test SafeBluesCampusSimulation.hour(150) == 6  # 150 -> 5am -- 6am

    @test SafeBluesCampusSimulation.day(1) == 1    # 1 -> Monday
    @test SafeBluesCampusSimulation.day(25) == 2   # 25 -> Tuesday
    @test SafeBluesCampusSimulation.day(50) == 3   # 50 -> Wednesday
    @test SafeBluesCampusSimulation.day(100) == 5  # 100 -> Friday
    @test SafeBluesCampusSimulation.day(150) == 7  # 150 -> Sunday

    @test SafeBluesCampusSimulation.is_weekend(25) == false   # 25 -> Tuesday
    @test SafeBluesCampusSimulation.is_weekend(50) == false   # 50 -> Wednesday
    @test SafeBluesCampusSimulation.is_weekend(100) == false  # 100 -> Friday
    @test SafeBluesCampusSimulation.is_weekend(125) == true   # 125 -> Saturday
    @test SafeBluesCampusSimulation.is_weekend(150) == true   # 150 -> Sunday
end
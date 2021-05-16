# Tests for time-related functions defined in `src/core/simulate.jl`.

@testset "Time" begin
    @test hour(1) == 1    # 1 -> 12am -- 1am
    @test hour(25) == 1   # 25 -> 12am -- 1am
    @test hour(50) == 2   # 50 -> 1am -- 2am
    @test hour(100) == 4  # 100 -> 3am -- 4am
    @test hour(150) == 6  # 150 -> 5am -- 6am

    @test day(1) == 1    # 1 -> Monday
    @test day(25) == 2   # 25 -> Tuesday
    @test day(50) == 3   # 50 -> Wednesday
    @test day(100) == 5  # 100 -> Friday
    @test day(150) == 7  # 150 -> Sunday

    @test is_weekend(25) == false   # 25 -> Tuesday
    @test is_weekend(50) == false   # 50 -> Wednesday
    @test is_weekend(100) == false  # 100 -> Friday
    @test is_weekend(125) == true   # 125 -> Saturday
    @test is_weekend(150) == true   # 150 -> Sunday
end
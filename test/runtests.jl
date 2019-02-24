using ModernRobotics
using Test

@testset "ModernRobotics.jl" begin
    @testset "basic helper functions" begin
        @test NearZero(-1e-7)
        @test Normalize([1, 2, 3]) == [0.2672612419124244,
                                       0.5345224838248488,
                                       0.8017837257372732]
    end
    @testset "chapter 3: rigid-body motions" begin
        @test RotInv([0 0 1; 1 0 0; 0 1 0]) == [0 1 0; 0 0 1; 1 0 0]
        @test VecToso3([1 2 3]) == [0 -3 2; 3 0 -1; -2 1 0]
        @test so3ToVec([0 -3 2; 3 0 -1; -2 1 0]) == [1, 2, 3]
        @testset "AxisAng3" begin
            omghat, θ = AxisAng3([1, 2, 3])
            @test omghat == [0.2672612419124244,
                             0.5345224838248488,
                             0.8017837257372732]
            @test θ == 3.7416573867739413
        end
    end
    @testset "chapter 4: forward kinematics" begin
    end
    @testset "chapter 5: velocity kinematics and statics" begin
    end
    @testset "chapter 6: inverse kinematics" begin
    end
    @testset "chapter 8: dynamics of open chains" begin
    end
    @testset "chapter 9: trajectory generation" begin
    end
    @testset "chapter 11: robot control" begin
    end
end

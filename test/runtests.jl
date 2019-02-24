using ModernRoboticsBook
using Test

@testset "ModernRoboticsBook.jl" begin
    @testset "basic helper functions" begin
        @test NearZero(-1e-7)
        @test Normalize([1, 2, 3]) == [0.2672612419124244,
                                       0.5345224838248488,
                                       0.8017837257372732]
    end
    @testset "chapter 3: rigid-body motions" begin
        @test RotInv([0 0 1; 1 0 0; 0 1 0]) == [0 1 0; 0 0 1; 1 0 0]
        @test VecToso3([1, 2, 3]) == [0 -3 2; 3 0 -1; -2 1 0]
        @test so3ToVec([0 -3 2; 3 0 -1; -2 1 0]) == [1, 2, 3]
        @test AxisAng3([1, 2, 3]) == ([0.2672612419124244,
                                       0.5345224838248488,
                                       0.8017837257372732], 3.7416573867739413)
        @test isapprox(MatrixExp3([ 0 -3  2;
                                    3  0 -1;
                                   -2  1  0]),
                       [-0.694921  0.713521 0.0892929
                        -0.192007 -0.303785 0.933192 
                         0.692978  0.63135  0.348107 ]; rtol=1e-6)
        @test isapprox(MatrixLog3([0 0 1;
                                   1 0 0;
                                   0 1 0]),
                       [ 0.0    -1.2092  1.2092
                         1.2092  0.0    -1.2092
                        -1.2092  1.2092  0.0   ]; rtol=1e-6)
        @test RpToTrans([1 0  0;
                         0 0 -1;
                         0 1  0], [1, 2, 5]) == [1 0  0 1
                                                 0 0 -1 2
                                                 0 1  0 5
                                                 0 0  0 1]
        @test TransToRp([1 0  0 0;
                         0 0 -1 0;
                         0 1  0 3;
                         0 0  0 1]) == ([1 0  0;
                                         0 0 -1;
                                         0 1  0], [0, 0, 3])
        @test TransInv([1 0  0 0;
                        0 0 -1 0;
                        0 1  0 3;
                        0 0  0 1]) == [1  0 0  0
                                       0  0 1 -3
                                       0 -1 0  0
                                       0  0 0  1]
        @test VecTose3([1, 2, 3, 4, 5, 6]) == [ 0 -3  2  4
                                                3  0 -1  5
                                               -2  1  0  6
                                                0  0  0  0]
        @test se3ToVec([ 0 -3  2  4;
                         3  0 -1  5;
                        -2  1  0  6;
                         0  0  0  0]) == [1, 2, 3, 4, 5, 6]
        @test Adjoint([1  0  0  0;
                       0  0 -1  0;
                       0  1  0  3;
                       0  0  0  1]) == [1  0  0  0  0  0
                                        0  0 -1  0  0  0
                                        0  1  0  0  0  0
                                        0  0  3  1  0  0
                                        3  0  0  0  0 -1
                                        0  0  0  0  1  0]
        @test ScrewToAxis([3; 0; 0], [0; 0; 1], 2) == [0, 0, 1, 0, -3, 2]
        @test AxisAng6([1, 0, 0, 1, 2, 3]) == ([1, 0, 0, 1, 2, 3], 1)
        @test MatrixExp6([0    0    0    0   ;
                          0    0   -π/2  3π/4;
                          0    π/2  0    3π/4;
                          0    0    0    0   ]) ≈ [1  0  0  0;
                                                   0  0 -1  0;
                                                   0  1  0  3;
                                                   0  0  0  1]
        @test MatrixLog6([1  0  0  0;
                          0  0 -1  0;
                          0  1  0  3;
                          0  0  0  1]) ≈ [0    0    0    0   ;
                                          0    0   -π/2  3π/4;
                                          0    π/2  0    3π/4;
                                          0    0    0    0   ]
        @test isapprox(ProjectToSO3([ 0.675  0.150  0.720;
                                      0.370  0.771 -0.511;
                                     -0.630  0.619  0.472]),
                       [ 0.679011  0.148945   0.718859
                         0.373207  0.773196  -0.512723
                        -0.632187  0.616428   0.469421]; rtol=1e-6)
        @test isapprox(ProjectToSE3([ 0.675  0.150  0.720  1.2;
                                      0.370  0.771 -0.511  5.4;
                                     -0.630  0.619  0.472  3.6;
                                      0.003  0.002  0.010  0.9]),
                                    [ 0.679011  0.148945   0.718859  1.2
                                      0.373207  0.773196  -0.512723  5.4
                                     -0.632187  0.616428   0.469421  3.6
                                      0.0       0.0        0.0       1.0]; rtol=1e-6)
        @test DistanceToSO3([1.0  0.0  0.0 ;
                             0.0  0.1 -0.95;
                             0.0  1.0  0.1 ]) == 0.08835298523536149
        @test DistanceToSE3([1.0  0.0  0.0   1.2 ;
                             0.0  0.1 -0.95  1.5 ;
                             0.0  1.0  0.1  -0.9 ;
                             0.0  0.0  0.1   0.98]) == 0.13493053768513638
        @test !TestIfSO3([1.0  0.0  0.0 ;
                          0.0  0.1 -0.95;
                          0.0  1.0  0.1 ])
        @test !TestIfSE3([1.0  0.0  0.0   1.2 ;
                          0.0  0.1 -0.95  1.5 ;
                          0.0  1.0  0.1  -0.9 ;
                          0.0  0.0  0.1   0.98])
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

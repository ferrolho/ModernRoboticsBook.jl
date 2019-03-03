using ModernRoboticsBook
using Test

import LinearAlgebra
const linalg = LinearAlgebra

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
        @test MatrixExp3(zeros(3, 3)) == [1 0 0; 0 1 0; 0 0 1]

        @test isapprox(MatrixLog3([0 0 1;
                                   1 0 0;
                                   0 1 0]),
                       [ 0.0    -1.2092  1.2092
                         1.2092  0.0    -1.2092
                        -1.2092  1.2092  0.0   ]; rtol=1e-6)
        @test MatrixLog3(zeros(3, 3)) == zeros(3, 3)

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
        M = [-1  0  0  0 ;
              0  1  0  6 ;
              0  0 -1  2 ;
              0  0  0  1 ]

        Blist = [ 0  0 -1  2  0  0   ;
                  0  0  0  0  1  0   ;
                  0  0  1  0  0  0.1 ]'

        Slist = [ 0  0  1  4  0  0   ;
                  0  0  0  0  1  0   ;
                  0  0 -1 -6  0 -0.1 ]'

        thetalist = [ π/2, 3, π ]

        @test isapprox(FKinBody(M, Blist, thetalist),
                       [-1.14424e-17  1.0           0.0  -5.0     ;
                         1.0          1.14424e-17   0.0   4.0     ;
                         0.0          0.0          -1.0   1.68584 ;
                         0.0          0.0           0.0   1.0     ]; rtol=1e-6)
        @test isapprox(FKinSpace(M, Slist, thetalist),
                       [-1.14424e-17  1.0           0.0  -5.0     ;
                         1.0          1.14424e-17   0.0   4.0     ;
                         0.0          0.0          -1.0   1.68584 ;
                         0.0          0.0           0.0   1.0     ]; rtol=1e-6)
    end
    @testset "chapter 5: velocity kinematics and statics" begin
        Blist = [0 0 1   0 0.2 0.2;
                 1 0 0   2   0   3;
                 0 1 0   0   2   1;
                 1 0 0 0.2 0.3 0.4]'

        Slist = [0 0 1   0 0.2 0.2;
                 1 0 0   2   0   3;
                 0 1 0   0   2   1;
                 1 0 0 0.2 0.3 0.4]'

        thetalist = [0.2, 1.1, 0.1, 1.2]

        @test isapprox(JacobianBody(Blist, thetalist),
                       [ -0.0452841  0.995004    0.0       1.0 ;
                          0.743593   0.0930486   0.362358  0.0 ;
                         -0.667097   0.0361754  -0.932039  0.0 ;
                          2.32586    1.66809     0.564108  0.2 ;
                         -1.44321    2.94561     1.43307   0.3 ;
                         -2.0664     1.82882    -1.58869   0.4 ]; rtol=1e-5)
        @test isapprox(JacobianSpace(Slist, thetalist),
                       [0.0  0.980067  -0.0901156   0.957494 ;
                        0.0  0.198669   0.444554    0.284876 ;
                        1.0  0.0        0.891207   -0.0452841;
                        0.0  1.95219   -2.21635    -0.511615 ;
                        0.2  0.436541  -2.43713     2.77536  ;
                        0.2  2.96027    3.23573     2.22512  ]; rtol=1e-5)
    end
    @testset "chapter 6: inverse kinematics" begin
        Blist = [ 0  0 -1  2  0  0   ;
                  0  0  0  0  1  0   ;
                  0  0  1  0  0  0.1 ]'

        Slist = [ 0  0  1  4  0  0   ;
                  0  0  0  0  1  0   ;
                  0  0 -1 -6  0 -0.1 ]'

        M = [ -1  0  0  0 ;
               0  1  0  6 ;
               0  0 -1  2 ;
               0  0  0  1 ]

        T = [ 0  1  0     -5 ;
              1  0  0      4 ;
              0  0 -1 1.6858 ;
              0  0  0      1 ]

        thetalist0 = [1.5, 2.5, 3]

        eomg, ev = 0.01, 0.001

        thetalist, success = IKinBody(Blist, M, T, thetalist0, eomg, ev)
        @test isapprox(thetalist, [1.57074, 2.99967, 3.14154]; rtol=1e-5)
        @test success

        thetalist, success = IKinSpace(Slist, M, T, thetalist0, eomg, ev)
        @test isapprox(thetalist, [1.57074, 2.99966, 3.14153]; rtol=1e-5)
        @test success
    end
    @testset "chapter 8: dynamics of open chains" begin
        @test ad([1, 2, 3, 4, 5, 6]) == [ 0.0  -3.0   2.0   0.0   0.0   0.0
                                          3.0   0.0  -1.0   0.0   0.0   0.0
                                         -2.0   1.0   0.0   0.0   0.0   0.0
                                          0.0  -6.0   5.0   0.0  -3.0   2.0
                                          6.0   0.0  -4.0   3.0   0.0  -1.0
                                         -5.0   4.0   0.0  -2.0   1.0   0.0 ]

        thetalist = [0.1, 0.1, 0.1]
        dthetalist = [0.1, 0.2, 0.3]
        ddthetalist = [2, 1.5, 1]
        g = [0, 0, -9.8]
        Ftip = [1, 1, 1, 1, 1, 1]
        M01 = [1  0  0         0 ;
               0  1  0         0 ;
               0  0  1  0.089159 ;
               0  0  0         1 ]
        M12 = [0  0  1     0.28 ;
               0  1  0  0.13585 ;
              -1  0  0        0 ;
               0  0  0        1 ]
        M23 = [1  0  0        0 ;
               0  1  0  -0.1197 ;
               0  0  1    0.395 ;
               0  0  0        1 ]
        M34 = [1  0  0        0 ;
               0  1  0        0 ;
               0  0  1  0.14225 ;
               0  0  0        1 ]
        G1 = linalg.Diagonal([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
        G2 = linalg.Diagonal([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
        G3 = linalg.Diagonal([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])
        Glist = [G1, G2, G3]
        Mlist = [M01, M12, M23, M34]
        Slist = [ 1  0  1      0  1      0 ;
                  0  1  0 -0.089  0      0 ;
                  0  1  0 -0.089  0  0.425 ]'

        taulist_actual = InverseDynamics(thetalist, dthetalist, ddthetalist, g, Ftip, Mlist, Glist, Slist)

        @test taulist_actual ≈ [74.69616155287451, -33.06766015851458, -3.230573137901424]

        @test isapprox(MassMatrix(thetalist, Mlist, Glist, Slist),
                       [ 22.5433      -0.307147  -0.00718426;
                         -0.307147     1.96851    0.432157  ;
                         -0.00718426   0.432157   0.191631  ]; rtol=1e-5)

        @test VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist) ≈ [  0.26453118054501235 ;
                                                                                 -0.0550515682891655  ;
                                                                                 -0.006891320068248911]

        @test GravityForces(thetalist, g, Mlist, Glist, Slist) == [  28.40331261821983  ;
                                                                    -37.64094817177068  ;
                                                                     -5.4415891999683605]

        @test EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist) ≈ [ 1.4095460782639782;
                                                                          1.8577149723180628;
                                                                          1.392409          ]

        taulist = [0.5, 0.6, 0.7]
        @test ForwardDynamics(thetalist, dthetalist, taulist, g, Ftip, Mlist, Glist, Slist) ≈ [ -0.9739290670855626;
                                                                                                25.584667840340558 ;
                                                                                               -32.91499212478149  ]

        @test hcat(EulerStep([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], [2, 1.5, 1], 0.1)...) ≈ hcat([0.11, 0.12, 0.13], [0.3, 0.35, 0.4])

        @testset "chapter 8: dynamics of open chains" begin
            thetastart = [0, 0, 0]
            thetaend = [π/2, π/2, π/2]
            Tf = 3
            N = 10
            method = 5

            traj = JointTrajectory(thetastart, thetaend, Tf, N, method)

            thetamat = copy(traj)
            dthetamat = zeros(N, 3)
            ddthetamat = zeros(N, 3)

            dt = Tf / (N - 1)

            for i = 1:size(traj, 1) - 1
                dthetamat[i + 1, :] = (thetamat[i + 1, :] - thetamat[i, :]) / dt
                ddthetamat[i + 1, :] = (dthetamat[i + 1, :] - dthetamat[i, :]) / dt
            end

            Ftipmat = ones(N, 6)

            taumat_actual = InverseDynamicsTrajectory(thetamat, dthetamat, ddthetamat, g, Ftipmat, Mlist, Glist, Slist)

            taumat_expected = [ 1.32297079e+01 -3.62621080e+01 -4.18134100e+00 ;
                                1.96974217e+01 -3.59188701e+01 -4.07919728e+00 ;
                                5.11979532e+01 -3.44705050e+01 -3.59488765e+00 ;
                                9.41368122e+01 -3.14099606e+01 -2.41622731e+00 ;
                                1.25417579e+02 -2.45283212e+01  5.76295281e-02 ;
                                1.24948454e+02 -1.85038921e+01  1.94898550e+00 ;
                                1.01525941e+02 -1.88375820e+01  2.07369432e+00 ;
                                7.68134579e+01 -2.23610568e+01  1.96172027e+00 ;
                                6.44365965e+01 -2.46704062e+01  2.07890001e+00 ;
                                6.72354909e+01 -2.47008371e+01  2.19474783e+00 ]

            @test taumat_actual ≈ taumat_expected
        end
    end
    @testset "chapter 9: trajectory generation" begin
        @test CubicTimeScaling(2, 0.6) == 0.21600000000000003
        @test QuinticTimeScaling(2, 0.6) == 0.16308
    end
    @testset "chapter 11: robot control" begin
        thetalist = [0.1, 0.1, 0.1]
        dthetalist = [0.1, 0.2, 0.3]

        # Initialise robot description (Example with 3 links)
        g = [0, 0, -9.8]

        M01 = [1  0  0         0 ;
               0  1  0         0 ;
               0  0  1  0.089159 ;
               0  0  0         1 ]
        M12 = [0  0  1     0.28 ;
               0  1  0  0.13585 ;
              -1  0  0        0 ;
               0  0  0        1 ]
        M23 = [1  0  0        0 ;
               0  1  0  -0.1197 ;
               0  0  1    0.395 ;
               0  0  0        1 ]
        M34 = [1  0  0        0 ;
               0  1  0        0 ;
               0  0  1  0.14225 ;
               0  0  0        1 ]

        G1 = linalg.Diagonal([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
        G2 = linalg.Diagonal([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
        G3 = linalg.Diagonal([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])

        Mlist = [M01, M12, M23, M34]
        Glist = [G1, G2, G3]

        Slist = [ 1  0  1      0  1      0 ;
                0  1  0 -0.089  0      0 ;
                0  1  0 -0.089  0  0.425 ]'

        # Create a trajectory to follow
        thetaend = [π / 2, π, 1.5 * π]

        Tf = 1
        dt = 0.05
        N = Int(Tf / dt)
        method = 5

        traj = JointTrajectory(thetalist, thetaend, Tf, N, method)

        thetamatd = copy(traj)
        dthetamatd = zeros(N, 3)
        ddthetamatd = zeros(N, 3)
        dt = Tf / (N - 1)

        for i = 1:size(traj, 1) - 1
            dthetamatd[i + 1, :] = (thetamatd[i + 1, :] - thetamatd[i, :]) / dt
            ddthetamatd[i + 1, :] = (dthetamatd[i + 1, :] - dthetamatd[i, :]) / dt
        end

        # Possibly wrong robot description (Example with 3 links)
        gtilde = [0.8, 0.2, -8.8]

        Mhat01 = [1 0 0   0;
                  0 1 0   0;
                  0 0 1 0.1;
                  0 0 0   1]
        Mhat12 = [0 0 1 0.3;
                  0 1 0 0.2;
                 -1 0 0   0;
                  0 0 0   1]
        Mhat23 = [1 0 0    0;
                  0 1 0 -0.2;
                  0 0 1  0.4;
                  0 0 0    1]
        Mhat34 = [1 0 0   0;
                  0 1 0   0;
                  0 0 1 0.2;
                  0 0 0   1]

        Ghat1 = linalg.Diagonal([0.1, 0.1, 0.1, 4, 4, 4])
        Ghat2 = linalg.Diagonal([0.3, 0.3, 0.1, 9, 9, 9])
        Ghat3 = linalg.Diagonal([0.1, 0.1, 0.1, 3, 3, 3])

        Gtildelist = [Ghat1, Ghat2, Ghat3]
        Mtildelist = [Mhat01, Mhat12, Mhat23, Mhat34]

        # Other required arguments
        Ftipmat = ones(size(traj, 1), 6)

        Kp = 20
        Ki = 10
        Kd = 18
        intRes = 8

        taumat_expected = [ -14.2640765   -54.06797429  -11.265448   ;
                             71.7014572   -17.58330542    3.86417108 ;
                            208.80692807    6.94442209    8.4352746  ;
                            269.9223766    14.44412677   11.24081382 ;
                            316.48343344    6.4020598    10.60970699 ;
                            327.82241593   -3.98984379   14.31752441 ;
                            248.33306921  -16.39336633   21.61795095 ;
                             93.7564835   -28.5575642    28.0092122  ;
                             13.12918592  -44.38407547   19.04258057 ;
                             56.35246455   -8.56189073    1.69770764 ;
                             32.68030349   39.77791901   -5.94800597 ;
                            -49.85502041   37.95496258  -15.10367806 ;
                           -104.48630504   24.8129766   -16.25667052 ;
                           -123.14920836   -3.62727714  -14.4680164  ;
                            -84.15220471  -25.40152665  -12.85439272 ;
                            -50.09890916  -33.73575763  -12.38441089 ;
                            -30.41466046  -34.03362524  -11.30974711 ;
                            -13.90701987  -26.40700004   -9.50026445 ;
                              7.93416317  -12.95474009   -6.58132646 ;
                             44.38627151    4.53534718   -2.32611269 ]

        thetamat_expected = [ 0.1028237  0.10738308 0.07715206 ;
                              0.10475535 0.11170712 0.05271794 ;
                              0.11930448 0.13796752 0.12315113 ;
                              0.1582068  0.21615938 0.27654789 ;
                              0.22464764 0.35314707 0.51117109 ;
                              0.31872944 0.54551689 0.81897456 ;
                              0.4377715  0.78011155 1.21554584 ;
                              0.57501633 1.03974866 1.71548388 ;
                              0.72137128 1.30551854 2.2916328  ;
                              0.87221934 1.56564069 2.84693996 ;
                              1.0257972  1.841874   3.32092634 ;
                              1.17352729 2.14460368 3.72255083 ;
                              1.30459739 2.44287634 4.03953656 ;
                              1.41155584 2.71108778 4.28079221 ;
                              1.49166143 2.92224332 4.45724092 ;
                              1.54707169 3.06503575 4.57357234 ;
                              1.579692   3.14295613 4.62867261 ;
                              1.59184568 3.1665363  4.62924012 ;
                              1.58827951 3.15316968 4.59115503 ;
                              1.57746305 3.12631855 4.54394792 ]

        taumat_actual, thetamat_actual = SimulateControl(thetalist, dthetalist, g, Ftipmat, Mlist, Glist,
                                                         Slist, thetamatd, dthetamatd, ddthetamatd, gtilde,
                                                         Mtildelist, Gtildelist, Kp, Ki, Kd, dt, intRes)

        @test taumat_actual ≈ taumat_expected
        @test thetamat_actual ≈ thetamat_expected
    end
end

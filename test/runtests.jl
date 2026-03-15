using Aqua
using JSON
using ModernRoboticsBook
using StaticArrays
using Test

import LinearAlgebra as LA

# Re-use the internal JSON matrix helper for parsing reference data
_json_matrix(rows) = reduce(vcat, [Float64.(r)' for r in rows])

Aqua.test_all(ModernRoboticsBook)

@testset "ModernRoboticsBook.jl" begin
    @testset "basic helper functions" begin
        @test near_zero(-1e-7)
    end
    @testset "chapter 3: rigid-body motions" begin
        @test vec_to_so3([1, 2, 3]) == [0 -3 2; 3 0 -1; -2 1 0]
        @test so3_to_vec([0 -3 2; 3 0 -1; -2 1 0]) == [1, 2, 3]
        @test axis_angle3([1, 2, 3]) == (
            [0.2672612419124244, 0.5345224838248488, 0.8017837257372732],
            3.7416573867739413,
        )

        @test isapprox(
            matrix_exp3([
                0 -3 2
                3 0 -1
                -2 1 0
            ]),
            [
                -0.694921 0.713521 0.0892929
                -0.192007 -0.303785 0.933192
                0.692978 0.63135 0.348107
            ];
            rtol = 1e-6,
        )
        @test matrix_exp3(zeros(3, 3)) == [1 0 0; 0 1 0; 0 0 1]

        @test isapprox(
            matrix_log3([
                0 0 1
                1 0 0
                0 1 0
            ]),
            [
                0.0 -1.2092 1.2092
                1.2092 0.0 -1.2092
                -1.2092 1.2092 0.0
            ];
            rtol = 1e-6,
        )
        @test matrix_log3(zeros(3, 3)) == zeros(3, 3)
        @test matrix_log3([-3 0 0; 0 1 0; 0 0 1]) == [0 -π 0; π 0 0; 0 0 0]
        @test matrix_log3([-2 0 0; 0 1 0; 0 0 -1]) == [0 0 π; 0 0 0; -π 0 0]
        @test matrix_log3([1 0 0; 0 -1 0; 0 0 -1]) == [0 0 0; 0 0 -π; 0 π 0]

        @test rotation_position_to_transform(
            [
                1 0 0
                0 0 -1
                0 1 0
            ],
            [1, 2, 5],
        ) == [
            1 0 0 1
            0 0 -1 2
            0 1 0 5
            0 0 0 1
        ]
        @test transform_to_rotation_position([
            1 0 0 0
            0 0 -1 0
            0 1 0 3
            0 0 0 1
        ]) == (
            [
                1 0 0
                0 0 -1
                0 1 0
            ],
            [0, 0, 3],
        )
        @test transform_inv([
            1 0 0 0
            0 0 -1 0
            0 1 0 3
            0 0 0 1
        ]) == [
            1 0 0 0
            0 0 1 -3
            0 -1 0 0
            0 0 0 1
        ]
        @test vec_to_se3([1, 2, 3, 4, 5, 6]) == [
            0 -3 2 4
            3 0 -1 5
            -2 1 0 6
            0 0 0 0
        ]
        @test se3_to_vec([
            0 -3 2 4
            3 0 -1 5
            -2 1 0 6
            0 0 0 0
        ]) == [1, 2, 3, 4, 5, 6]
        @test adjoint_representation([
            1 0 0 0
            0 0 -1 0
            0 1 0 3
            0 0 0 1
        ]) == [
            1 0 0 0 0 0
            0 0 -1 0 0 0
            0 1 0 0 0 0
            0 0 3 1 0 0
            3 0 0 0 0 -1
            0 0 0 0 1 0
        ]
        @test screw_to_axis([3; 0; 0], [0; 0; 1], 2) == [0, 0, 1, 0, -3, 2]

        @test axis_angle6([1, 0, 0, 1, 2, 3]) == ([1, 0, 0, 1, 2, 3], 1)
        @test axis_angle6([0, 0, 0, 0, 0, 4]) == ([0, 0, 0, 0, 0, 1], 4)

        @test matrix_exp6([
            0 0 0 0
            0 0 -π/2 3π/4
            0 π/2 0 3π/4
            0 0 0 0
        ]) ≈ [
            1 0 0 0
            0 0 -1 0
            0 1 0 3
            0 0 0 1
        ]
        @test matrix_log6([
            1 0 0 0
            0 0 -1 0
            0 1 0 3
            0 0 0 1
        ]) ≈ [
            0 0 0 0
            0 0 -π/2 3π/4
            0 π/2 0 3π/4
            0 0 0 0
        ]
        @test matrix_log6(Array(LA.Diagonal(ones(4)))) == zeros(4, 4)

        @test isapprox(
            project_to_so3([
                0.675 0.150 0.720
                0.370 0.771 -0.511
                -0.630 0.619 0.472
            ]),
            [
                0.679011 0.148945 0.718859
                0.373207 0.773196 -0.512723
                -0.632187 0.616428 0.469421
            ];
            rtol = 1e-6,
        )
        @test project_to_so3([-1 0 0; 0 -1 0; 0 0 -1]) == [-1 0 0; 0 -1 0; 0 0 1]

        @test isapprox(
            project_to_se3(
                [
                    0.675 0.150 0.720 1.2
                    0.370 0.771 -0.511 5.4
                    -0.630 0.619 0.472 3.6
                    0.003 0.002 0.010 0.9
                ],
            ),
            [
                0.679011 0.148945 0.718859 1.2
                0.373207 0.773196 -0.512723 5.4
                -0.632187 0.616428 0.469421 3.6
                0.0 0.0 0.0 1.0
            ];
            rtol = 1e-6,
        )

        @test distance_to_so3([
            1.0 0.0 0.0
            0.0 0.1 -0.95
            0.0 1.0 0.1
        ]) == 0.08835298523536149

        @test distance_to_se3(
            [
                1.0 0.0 0.0 1.2
                0.0 0.1 -0.95 1.5
                0.0 1.0 0.1 -0.9
                0.0 0.0 0.1 0.98
            ],
        ) == 0.13493053768513638
        @test distance_to_se3(Array(reshape(1:9, (3, 3)))) >= 1e+9

        @test !is_so3([
            1.0 0.0 0.0
            0.0 0.1 -0.95
            0.0 1.0 0.1
        ])
        @test !is_se3([
            1.0 0.0 0.0 1.2
            0.0 0.1 -0.95 1.5
            0.0 1.0 0.1 -0.9
            0.0 0.0 0.1 0.98
        ])
    end
    @testset "chapter 4: forward kinematics" begin
        home_ee_pose = [
            -1 0 0 0
            0 1 0 6
            0 0 -1 2
            0 0 0 1
        ]

        body_screw_axes = [
            0 0 -1 2 0 0
            0 0 0 0 1 0
            0 0 1 0 0 0.1
        ]'

        screw_axes = [
            0 0 1 4 0 0
            0 0 0 0 1 0
            0 0 -1 -6 0 -0.1
        ]'

        joint_positions = [π / 2, 3, π]

        @test isapprox(
            forward_kinematics_body(home_ee_pose, body_screw_axes, joint_positions),
            [
                -1.14424e-17 1.0 0.0 -5.0
                1.0 1.14424e-17 0.0 4.0
                0.0 0.0 -1.0 1.68584
                0.0 0.0 0.0 1.0
            ];
            rtol = 1e-6,
        )
        @test isapprox(
            forward_kinematics_space(home_ee_pose, screw_axes, joint_positions),
            [
                -1.14424e-17 1.0 0.0 -5.0
                1.0 1.14424e-17 0.0 4.0
                0.0 0.0 -1.0 1.68584
                0.0 0.0 0.0 1.0
            ];
            rtol = 1e-6,
        )
    end
    @testset "chapter 5: velocity kinematics and statics" begin
        body_screw_axes = [
            0 0 1 0 0.2 0.2
            1 0 0 2 0 3
            0 1 0 0 2 1
            1 0 0 0.2 0.3 0.4
        ]'

        screw_axes = [
            0 0 1 0 0.2 0.2
            1 0 0 2 0 3
            0 1 0 0 2 1
            1 0 0 0.2 0.3 0.4
        ]'

        joint_positions = [0.2, 1.1, 0.1, 1.2]

        @test isapprox(
            jacobian_body(body_screw_axes, joint_positions),
            [
                -0.0452841 0.995004 0.0 1.0
                0.743593 0.0930486 0.362358 0.0
                -0.667097 0.0361754 -0.932039 0.0
                2.32586 1.66809 0.564108 0.2
                -1.44321 2.94561 1.43307 0.3
                -2.0664 1.82882 -1.58869 0.4
            ];
            rtol = 1e-5,
        )
        @test isapprox(
            jacobian_space(screw_axes, joint_positions),
            [
                0.0 0.980067 -0.0901156 0.957494
                0.0 0.198669 0.444554 0.284876
                1.0 0.0 0.891207 -0.0452841
                0.0 1.95219 -2.21635 -0.511615
                0.2 0.436541 -2.43713 2.77536
                0.2 2.96027 3.23573 2.22512
            ];
            rtol = 1e-5,
        )
    end
    @testset "chapter 6: inverse kinematics" begin
        body_screw_axes = [
            0 0 -1 2 0 0
            0 0 0 0 1 0
            0 0 1 0 0 0.1
        ]'

        screw_axes = [
            0 0 1 4 0 0
            0 0 0 0 1 0
            0 0 -1 -6 0 -0.1
        ]'

        home_ee_pose = [
            -1 0 0 0
            0 1 0 6
            0 0 -1 2
            0 0 0 1
        ]

        target_config = [
            0 1 0 -5
            1 0 0 4
            0 0 -1 1.6858
            0 0 0 1
        ]

        initial_guess = [1.5, 2.5, 3]

        angular_tolerance, linear_tolerance = 0.01, 0.001

        joint_positions, success = inverse_kinematics_body(
            body_screw_axes,
            home_ee_pose,
            target_config,
            initial_guess,
            angular_tolerance,
            linear_tolerance,
        )
        @test isapprox(joint_positions, [1.57074, 2.99967, 3.14154]; rtol = 1e-5)
        @test success

        joint_positions, success = inverse_kinematics_space(
            screw_axes,
            home_ee_pose,
            target_config,
            initial_guess,
            angular_tolerance,
            linear_tolerance,
        )
        @test isapprox(joint_positions, [1.57074, 2.99966, 3.14153]; rtol = 1e-5)
        @test success
    end
    @testset "chapter 8: dynamics of open chains" begin
        @test ad([1, 2, 3, 4, 5, 6]) == [
            0.0 -3.0 2.0 0.0 0.0 0.0
            3.0 0.0 -1.0 0.0 0.0 0.0
            -2.0 1.0 0.0 0.0 0.0 0.0
            0.0 -6.0 5.0 0.0 -3.0 2.0
            6.0 0.0 -4.0 3.0 0.0 -1.0
            -5.0 4.0 0.0 -2.0 1.0 0.0
        ]

        joint_positions = [0.1, 0.1, 0.1]
        joint_velocities = [0.1, 0.2, 0.3]
        joint_accelerations = [2, 1.5, 1]
        gravity = [0, 0, -9.8]
        tip_wrench = [1, 1, 1, 1, 1, 1]
        M01 = [
            1 0 0 0
            0 1 0 0
            0 0 1 0.089159
            0 0 0 1
        ]
        M12 = [
            0 0 1 0.28
            0 1 0 0.13585
            -1 0 0 0
            0 0 0 1
        ]
        M23 = [
            1 0 0 0
            0 1 0 -0.1197
            0 0 1 0.395
            0 0 0 1
        ]
        M34 = [
            1 0 0 0
            0 1 0 0
            0 0 1 0.14225
            0 0 0 1
        ]
        G1 = LA.Diagonal([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
        G2 = LA.Diagonal([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
        G3 = LA.Diagonal([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])
        spatial_inertias = [G1, G2, G3]
        link_frames = [M01, M12, M23, M34]
        screw_axes = [
            1 0 1 0 1 0
            0 1 0 -0.089 0 0
            0 1 0 -0.089 0 0.425
        ]'

        joint_torques_actual = inverse_dynamics_rnea(
            joint_positions,
            joint_velocities,
            joint_accelerations,
            gravity,
            tip_wrench,
            link_frames,
            spatial_inertias,
            screw_axes,
        )

        @test joint_torques_actual ≈
              [74.69616155287451, -33.06766015851458, -3.230573137901424]

        @test isapprox(
            mass_matrix_crba(joint_positions, link_frames, spatial_inertias, screw_axes),
            [
                22.5433 -0.307147 -0.00718426
                -0.307147 1.96851 0.432157
                -0.00718426 0.432157 0.191631
            ];
            rtol = 1e-5,
        )

        @test velocity_quadratic_forces(
            joint_positions,
            joint_velocities,
            link_frames,
            spatial_inertias,
            screw_axes,
        ) ≈ [
            0.26453118054501235
            -0.0550515682891655
            -0.006891320068248911
        ]

        @test gravity_forces(
            joint_positions,
            gravity,
            link_frames,
            spatial_inertias,
            screw_axes,
        ) ≈ [
            28.40331261821983
            -37.64094817177068
            -5.4415891999683605
        ]

        @test end_effector_forces(
            joint_positions,
            tip_wrench,
            link_frames,
            spatial_inertias,
            screw_axes,
        ) ≈ [
            1.4095460782639782
            1.8577149723180628
            1.392409
        ]

        joint_torques = [0.5, 0.6, 0.7]
        @test forward_dynamics_crba(
            joint_positions,
            joint_velocities,
            joint_torques,
            gravity,
            tip_wrench,
            link_frames,
            spatial_inertias,
            screw_axes,
        ) ≈ [
            -0.9739290670855626
            25.584667840340558
            -32.91499212478149
        ]

        @test hcat(euler_step([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], [2, 1.5, 1], 0.1)...) ≈
              hcat([0.11, 0.12, 0.13], [0.3, 0.35, 0.4])

        @testset "inverse dynamics trajectory" begin
            joint_position_start = [0, 0, 0]
            joint_position_end = [π / 2, π / 2, π / 2]
            total_time = 3
            N = 10
            method = 5

            traj = joint_trajectory(
                joint_position_start,
                joint_position_end,
                total_time,
                N,
                method,
            )

            joint_position_traj = copy(traj)
            joint_velocity_traj = zeros(N, 3)
            joint_acceleration_traj = zeros(N, 3)

            timestep = total_time / (N - 1)

            for i = 1:(size(traj, 1)-1)
                joint_velocity_traj[i+1, :] =
                    (joint_position_traj[i+1, :] - joint_position_traj[i, :]) / timestep
                joint_acceleration_traj[i+1, :] =
                    (joint_velocity_traj[i+1, :] - joint_velocity_traj[i, :]) / timestep
            end

            tip_wrench_traj = ones(N, 6)

            joint_torque_traj_actual = inverse_dynamics_trajectory(
                joint_position_traj,
                joint_velocity_traj,
                joint_acceleration_traj,
                gravity,
                tip_wrench_traj,
                link_frames,
                spatial_inertias,
                screw_axes,
            )

            joint_torque_traj_expected = [
                1.32297079e+01 -3.62621080e+01 -4.18134100e+00
                1.96974217e+01 -3.59188701e+01 -4.07919728e+00
                5.11979532e+01 -3.44705050e+01 -3.59488765e+00
                9.41368122e+01 -3.14099606e+01 -2.41622731e+00
                1.25417579e+02 -2.45283212e+01 5.76295281e-02
                1.24948454e+02 -1.85038921e+01 1.94898550e+00
                1.01525941e+02 -1.88375820e+01 2.07369432e+00
                7.68134579e+01 -2.23610568e+01 1.96172027e+00
                6.44365965e+01 -2.46704062e+01 2.07890001e+00
                6.72354909e+01 -2.47008371e+01 2.19474783e+00
            ]

            @test joint_torque_traj_actual ≈ joint_torque_traj_expected
        end
        @testset "forward dynamics trajectory" begin
            joint_torque_traj = [
                3.63 -6.58 -5.57
                3.74 -5.55 -5.5
                4.31 -0.68 -5.19
                5.18 5.63 -4.31
                5.85 8.17 -2.59
                5.78 2.79 -1.7
                4.99 -5.3 -1.19
                4.08 -9.41 0.07
                3.56 -10.1 0.97
                3.49 -9.41 1.23
            ]

            tip_wrench_traj = ones(size(joint_torque_traj, 1), 6)

            timestep = 0.1
            integration_resolution = 8

            joint_position_traj_expected = [
                0.1 0.1 0.1
                0.10643138 0.2625997 -0.22664947
                0.10197954 0.71581297 -1.22521632
                0.0801044 1.33930884 -2.28074132
                0.0282165 2.11957376 -3.07544297
                -0.07123855 2.87726666 -3.83289684
                -0.20136466 3.397858 -4.83821609
                -0.32380092 3.73338535 -5.98695747
                -0.41523262 3.85883317 -7.01130559
                -0.4638099 3.63178793 -7.63190052
            ]

            joint_velocity_traj_expected = [
                0.1 0.2 0.3
                0.01212502 3.42975773 -7.74792602
                -0.13052771 5.55997471 -11.22722784
                -0.35521041 7.11775879 -9.18173035
                -0.77358795 8.17307573 -7.05744594
                -1.2350231 6.35907497 -8.99784746
                -1.31426299 4.07685875 -11.18480509
                -1.06794821 2.49227786 -11.69748583
                -0.70264871 -0.55925705 -8.16067131
                -0.1455669 -4.57149985 -3.43135114
            ]

            joint_position_traj_actual, joint_velocity_traj_actual =
                forward_dynamics_trajectory(
                    joint_positions,
                    joint_velocities,
                    joint_torque_traj,
                    gravity,
                    tip_wrench_traj,
                    link_frames,
                    spatial_inertias,
                    screw_axes,
                    timestep,
                    integration_resolution,
                )

            @test joint_position_traj_actual ≈ joint_position_traj_expected
            @test joint_velocity_traj_actual ≈ joint_velocity_traj_expected
        end
    end
    @testset "chapter 9: trajectory generation" begin
        @test cubic_time_scaling(2, 0.6) ≈ 0.216
        @test quintic_time_scaling(2, 0.6) ≈ 0.16308
        @testset "joint trajectory" begin
            joint_position_start = [1, 0, 0, 1, 1, 0.2, 0, 1]
            joint_position_end = [1.2, 0.5, 0.6, 1.1, 2, 2, 0.9, 1]
            total_time = 4
            N = 6
            method = 3

            expected = [
                1 0 0 1 1 0.2 0 1
                1.0208 0.052 0.0624 1.0104 1.104 0.3872 0.0936 1
                1.0704 0.176 0.2112 1.0352 1.352 0.8336 0.3168 1
                1.1296 0.324 0.3888 1.0648 1.648 1.3664 0.5832 1
                1.1792 0.448 0.5376 1.0896 1.896 1.8128 0.8064 1
                1.2 0.5 0.6 1.1 2 2 0.9 1
            ]

            @test joint_trajectory(
                joint_position_start,
                joint_position_end,
                total_time,
                N,
                method,
            ) ≈ expected
        end
        @testset "screw trajectory" begin
            transform_start = [
                1 0 0 1
                0 1 0 0
                0 0 1 1
                0 0 0 1
            ]
            transform_end = [
                0 0 1 0.1
                1 0 0 0
                0 1 0 4.1
                0 0 0 1
            ]

            total_time = 5
            N = 4
            method = 3

            expected = [
                [
                    1 0 0 1
                    0 1 0 0
                    0 0 1 1
                    0 0 0 1
                ],
                [
                    0.904 -0.25 0.346 0.441
                    0.346 0.904 -0.25 0.529
                    -0.25 0.346 0.904 1.601
                    0 0 0 1
                ],
                [
                    0.346 -0.25 0.904 -0.117
                    0.904 0.346 -0.25 0.473
                    -0.25 0.904 0.346 3.274
                    0 0 0 1
                ],
                [
                    0 0 1 0.1
                    1 0 0 0
                    0 1 0 4.1
                    0 0 0 1
                ],
            ]

            actual = screw_trajectory(transform_start, transform_end, total_time, N, method)

            @test isapprox(actual, expected; rtol = 1e-3)
        end
        @testset "cartesian trajectory" begin
            transform_start = [
                1 0 0 1
                0 1 0 0
                0 0 1 1
                0 0 0 1
            ]
            transform_end = [
                0 0 1 0.1
                1 0 0 0
                0 1 0 4.1
                0 0 0 1
            ]

            total_time = 5
            N = 4
            method = 5


            expected = [
                [
                    1 0 0 1
                    0 1 0 0
                    0 0 1 1
                    0 0 0 1
                ],
                [
                    0.937 -0.214 0.277 0.811
                    0.277 0.937 -0.214 0
                    -0.214 0.277 0.937 1.651
                    0 0 0 1
                ],
                [
                    0.277 -0.214 0.937 0.289
                    0.937 0.277 -0.214 0
                    -0.214 0.937 0.277 3.449
                    0 0 0 1
                ],
                [
                    0 0 1 0.1
                    1 0 0 0
                    0 1 0 4.1
                    0 0 0 1
                ],
            ]

            actual =
                cartesian_trajectory(transform_start, transform_end, total_time, N, method)

            @test isapprox(actual, expected; rtol = 1e-3)
        end
    end
    @testset "chapter 11: robot control" begin
        joint_positions = [0.1, 0.1, 0.1]
        joint_velocities = [0.1, 0.2, 0.3]

        # Initialise robot description (Example with 3 links)
        gravity = [0, 0, -9.8]

        M01 = [
            1 0 0 0
            0 1 0 0
            0 0 1 0.089159
            0 0 0 1
        ]
        M12 = [
            0 0 1 0.28
            0 1 0 0.13585
            -1 0 0 0
            0 0 0 1
        ]
        M23 = [
            1 0 0 0
            0 1 0 -0.1197
            0 0 1 0.395
            0 0 0 1
        ]
        M34 = [
            1 0 0 0
            0 1 0 0
            0 0 1 0.14225
            0 0 0 1
        ]

        G1 = LA.Diagonal([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
        G2 = LA.Diagonal([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
        G3 = LA.Diagonal([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])

        link_frames = [M01, M12, M23, M34]
        spatial_inertias = [G1, G2, G3]

        screw_axes = [
            1 0 1 0 1 0
            0 1 0 -0.089 0 0
            0 1 0 -0.089 0 0.425
        ]'

        # Create a trajectory to follow
        joint_position_end = [π / 2, π, 1.5 * π]

        total_time = 1
        timestep = 0.05
        N = Int(total_time / timestep)
        method = 5

        traj = joint_trajectory(joint_positions, joint_position_end, total_time, N, method)

        desired_joint_position_traj = copy(traj)
        desired_joint_velocity_traj = zeros(N, 3)
        desired_joint_acceleration_traj = zeros(N, 3)
        timestep = total_time / (N - 1)

        for i = 1:(size(traj, 1)-1)
            desired_joint_velocity_traj[i+1, :] =
                (desired_joint_position_traj[i+1, :] - desired_joint_position_traj[i, :]) /
                timestep
            desired_joint_acceleration_traj[i+1, :] =
                (desired_joint_velocity_traj[i+1, :] - desired_joint_velocity_traj[i, :]) /
                timestep
        end

        # Possibly wrong robot description (Example with 3 links)
        estimated_gravity = [0.8, 0.2, -8.8]

        Mhat01 = [
            1 0 0 0
            0 1 0 0
            0 0 1 0.1
            0 0 0 1
        ]
        Mhat12 = [
            0 0 1 0.3
            0 1 0 0.2
            -1 0 0 0
            0 0 0 1
        ]
        Mhat23 = [
            1 0 0 0
            0 1 0 -0.2
            0 0 1 0.4
            0 0 0 1
        ]
        Mhat34 = [
            1 0 0 0
            0 1 0 0
            0 0 1 0.2
            0 0 0 1
        ]

        Ghat1 = LA.Diagonal([0.1, 0.1, 0.1, 4, 4, 4])
        Ghat2 = LA.Diagonal([0.3, 0.3, 0.1, 9, 9, 9])
        Ghat3 = LA.Diagonal([0.1, 0.1, 0.1, 3, 3, 3])

        estimated_spatial_inertias = [Ghat1, Ghat2, Ghat3]
        estimated_link_frames = [Mhat01, Mhat12, Mhat23, Mhat34]

        # Other required arguments
        tip_wrench_traj = ones(size(traj, 1), 6)

        Kp = 20
        Ki = 10
        Kd = 18
        integration_resolution = 8

        joint_torque_traj_expected = [
            -14.2640765 -54.06797429 -11.265448
            71.7014572 -17.58330542 3.86417108
            208.80692807 6.94442209 8.4352746
            269.9223766 14.44412677 11.24081382
            316.48343344 6.4020598 10.60970699
            327.82241593 -3.98984379 14.31752441
            248.33306921 -16.39336633 21.61795095
            93.7564835 -28.5575642 28.0092122
            13.12918592 -44.38407547 19.04258057
            56.35246455 -8.56189073 1.69770764
            32.68030349 39.77791901 -5.94800597
            -49.85502041 37.95496258 -15.10367806
            -104.48630504 24.8129766 -16.25667052
            -123.14920836 -3.62727714 -14.4680164
            -84.15220471 -25.40152665 -12.85439272
            -50.09890916 -33.73575763 -12.38441089
            -30.41466046 -34.03362524 -11.30974711
            -13.90701987 -26.40700004 -9.50026445
            7.93416317 -12.95474009 -6.58132646
            44.38627151 4.53534718 -2.32611269
        ]

        joint_position_traj_expected = [
            0.1028237 0.10738308 0.07715206
            0.10475535 0.11170712 0.05271794
            0.11930448 0.13796752 0.12315113
            0.1582068 0.21615938 0.27654789
            0.22464764 0.35314707 0.51117109
            0.31872944 0.54551689 0.81897456
            0.4377715 0.78011155 1.21554584
            0.57501633 1.03974866 1.71548388
            0.72137128 1.30551854 2.2916328
            0.87221934 1.56564069 2.84693996
            1.0257972 1.841874 3.32092634
            1.17352729 2.14460368 3.72255083
            1.30459739 2.44287634 4.03953656
            1.41155584 2.71108778 4.28079221
            1.49166143 2.92224332 4.45724092
            1.54707169 3.06503575 4.57357234
            1.579692 3.14295613 4.62867261
            1.59184568 3.1665363 4.62924012
            1.58827951 3.15316968 4.59115503
            1.57746305 3.12631855 4.54394792
        ]

        joint_torque_traj_actual, joint_position_traj_actual = simulate_control(
            joint_positions,
            joint_velocities,
            gravity,
            tip_wrench_traj,
            link_frames,
            spatial_inertias,
            screw_axes,
            desired_joint_position_traj,
            desired_joint_velocity_traj,
            desired_joint_acceleration_traj,
            estimated_gravity,
            estimated_link_frames,
            estimated_spatial_inertias,
            Kp,
            Ki,
            Kd,
            timestep,
            integration_resolution,
        )

        @test joint_torque_traj_actual ≈ joint_torque_traj_expected
        @test joint_position_traj_actual ≈ joint_position_traj_expected
    end

    @testset "Robot struct and load_robot" begin
        models_dir = joinpath(@__DIR__, "..", "robots", "models")

        @testset "load 3-DOF textbook model" begin
            robot = load_robot(joinpath(models_dir, "3dof_textbook.json"))
            @test robot.name == "3dof_textbook"
            @test robot.n_joints == 3
            @test size(robot.home_ee_pose) == (4, 4)
            @test size(robot.screw_axes_space) == (6, 3)
            @test size(robot.screw_axes_body) == (6, 3)
            @test length(robot.link_frames) == 4  # n+1
            @test length(robot.spatial_inertias) == 3
            @test robot.joint_types == [:revolute, :revolute, :revolute]
        end

        @testset "3-DOF Robot convenience wrappers match raw functions" begin
            robot = load_robot(joinpath(models_dir, "3dof_textbook.json"))

            q = [0.1, 0.1, 0.1]
            dq = [0.1, 0.2, 0.3]
            ddq = [2.0, 1.5, 1.0]
            g = [0.0, 0.0, -9.8]
            Ftip = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

            # Inverse dynamics via Robot wrapper must match raw function
            tau_robot =
                inverse_dynamics_rnea(robot, q, dq, ddq; gravity = g, tip_wrench = Ftip)
            tau_raw = inverse_dynamics_rnea(
                q,
                dq,
                ddq,
                g,
                Ftip,
                robot.link_frames,
                robot.spatial_inertias,
                robot.screw_axes_space,
            )
            @test tau_robot ≈ tau_raw

            # Mass matrix
            @test mass_matrix_crba(robot, q) ≈ mass_matrix_crba(
                q,
                robot.link_frames,
                robot.spatial_inertias,
                robot.screw_axes_space,
            )

            # Velocity quadratic forces
            @test velocity_quadratic_forces(robot, q, dq) ≈ velocity_quadratic_forces(
                q,
                dq,
                robot.link_frames,
                robot.spatial_inertias,
                robot.screw_axes_space,
            )

            # Gravity forces
            @test gravity_forces(robot, q; gravity = g) ≈ gravity_forces(
                q,
                g,
                robot.link_frames,
                robot.spatial_inertias,
                robot.screw_axes_space,
            )

            # End-effector forces
            @test end_effector_forces(robot, q, Ftip) ≈ end_effector_forces(
                q,
                Ftip,
                robot.link_frames,
                robot.spatial_inertias,
                robot.screw_axes_space,
            )

            # Forward dynamics
            joint_torques = [0.5, 0.6, 0.7]
            @test forward_dynamics_crba(
                robot,
                q,
                dq,
                joint_torques;
                gravity = g,
                tip_wrench = Ftip,
            ) ≈ forward_dynamics_crba(
                q,
                dq,
                joint_torques,
                g,
                Ftip,
                robot.link_frames,
                robot.spatial_inertias,
                robot.screw_axes_space,
            )
        end

        @testset "3-DOF dynamics match expected values" begin
            robot = load_robot(joinpath(models_dir, "3dof_textbook.json"))

            q = [0.1, 0.1, 0.1]
            dq = [0.1, 0.2, 0.3]
            ddq = [2.0, 1.5, 1.0]
            g = [0.0, 0.0, -9.8]
            Ftip = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

            @test inverse_dynamics_rnea(robot, q, dq, ddq; gravity = g, tip_wrench = Ftip) ≈
                  [74.69616155287451, -33.06766015851458, -3.230573137901424]

            @test isapprox(
                mass_matrix_crba(robot, q),
                [
                    22.5433 -0.307147 -0.00718426
                    -0.307147 1.96851 0.432157
                    -0.00718426 0.432157 0.191631
                ];
                rtol = 1e-5,
            )
        end

        @testset "textbook algorithms match optimized versions" begin
            robot = load_robot(joinpath(models_dir, "3dof_textbook.json"))

            q = [0.1, 0.1, 0.1]
            dq = [0.1, 0.2, 0.3]
            tau = [0.5, 0.6, 0.7]
            g = [0.0, 0.0, -9.8]
            Ftip = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

            @test mass_matrix_rnea(
                q,
                robot.link_frames,
                robot.spatial_inertias,
                robot.screw_axes_space,
            ) ≈ mass_matrix_crba(robot, q)

            @test forward_dynamics_rnea(
                q,
                dq,
                tau,
                g,
                Ftip,
                robot.link_frames,
                robot.spatial_inertias,
                robot.screw_axes_space,
            ) ≈ forward_dynamics_crba(robot, q, dq, tau; gravity = g, tip_wrench = Ftip)
        end

        @testset "load UR5 model" begin
            robot = load_robot(joinpath(models_dir, "ur5.json"))
            @test robot.name == "ur5"
            @test robot.n_joints == 6
            @test size(robot.screw_axes_space) == (6, 6)
            @test length(robot.link_frames) == 7  # n+1

            # FK at zero config should equal home config
            q0 = zeros(6)
            @test forward_kinematics_space(robot, q0) ≈ robot.home_ee_pose
        end

        @testset "load_robot by symbol" begin
            robot = load_robot(:ur5)
            @test robot.name == "ur5"
            @test robot.n_joints == 6
        end

        @testset "load_robot errors on missing model" begin
            @test_throws ErrorException load_robot(:nonexistent_robot)
        end

        @testset "zero allocations for core helpers" begin
            # Inputs as regular Matrix/Vector
            so3 = [0.0 -3.0 2.0; 3.0 0.0 -1.0; -2.0 1.0 0.0]
            R = [
                -0.694921 0.713521 0.0892929
                -0.192007 -0.303785 0.933192
                0.692978 0.63135 0.348107
            ]
            se3 = [
                0.0 -3.0 2.0 4.0
                3.0 0.0 -1.0 5.0
                -2.0 1.0 0.0 6.0
                0.0 0.0 0.0 0.0
            ]
            T = [
                1.0 0.0 0.0 0.0
                0.0 0.0 -1.0 0.0
                0.0 1.0 0.0 3.0
                0.0 0.0 0.0 1.0
            ]
            v3 = [1.0, 2.0, 3.0]
            v6 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]

            # Inputs as SMatrix/SVector
            sso3 = SMatrix{3,3}(so3)
            sR = SMatrix{3,3}(R)
            sse3 = SMatrix{4,4}(se3)
            sT = SMatrix{4,4}(T)
            sv3 = SVector{3}(v3)
            sv6 = SVector{6}(v6)

            # Test zero allocations with SMatrix/SVector inputs
            @test @allocated(vec_to_so3(sv3)) == 0
            @test @allocated(so3_to_vec(sso3)) == 0
            @test @allocated(matrix_exp3(sso3)) == 0
            @test @allocated(matrix_log3(sR)) == 0
            @test @allocated(vec_to_se3(sv6)) == 0
            @test @allocated(se3_to_vec(sse3)) == 0
            @test @allocated(transform_inv(sT)) == 0
            @test @allocated(adjoint_representation(sT)) == 0
            @test @allocated(matrix_exp6(sse3)) == 0
            @test @allocated(matrix_log6(sT)) == 0
            @test @allocated(ad(sv6)) == 0

            # Test zero allocations with regular Matrix/Vector inputs
            @test @allocated(vec_to_so3(v3)) == 0
            @test @allocated(so3_to_vec(so3)) == 0
            @test @allocated(vec_to_se3(v6)) == 0
            @test @allocated(se3_to_vec(se3)) == 0
            @test @allocated(ad(v6)) == 0
        end

        @testset "cross-validate UR5 against Pinocchio" begin
            # Pinocchio reference uses gravity [0, 0, -9.81] which matches our default
            robot = load_robot(:ur5)
            ref_path = joinpath(models_dir, "ur5_pinocchio_reference.json")
            ref = JSON.parsefile(ref_path)

            for (config_name, cfg) in ref["configurations"]
                @testset "$config_name" begin
                    q = Float64.(cfg["q"])
                    v = Float64.(cfg["v"])
                    a = Float64.(cfg["a"])

                    T_ref = _json_matrix(cfg["ee_pose"])
                    M_ref = _json_matrix(cfg["mass_matrix"])
                    g_ref = Float64.(cfg["gravity_forces"])
                    tau_ref = Float64.(cfg["inverse_dynamics"])
                    c_ref = Float64.(cfg["coriolis_forces"])

                    @test forward_kinematics_space(robot, q) ≈ T_ref atol = 1e-6
                    @test mass_matrix_crba(robot, q) ≈ M_ref atol = 1e-6
                    @test gravity_forces(robot, q) ≈ g_ref atol = 1e-6
                    @test inverse_dynamics_rnea(robot, q, v, a) ≈ tau_ref atol = 1e-6
                    @test velocity_quadratic_forces(robot, q, v) ≈ c_ref atol = 1e-6
                end
            end
        end
    end
end

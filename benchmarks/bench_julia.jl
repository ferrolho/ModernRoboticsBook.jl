#!/usr/bin/env julia
"""Benchmark ModernRoboticsBook.jl on UR5 for comparison with other libraries.

Run from the repo root:
    julia --project benchmarks/bench_julia.jl

Requires BenchmarkTools.jl (installed globally, not a project dependency):
    julia -e 'using Pkg; Pkg.add("BenchmarkTools")'
"""

using ModernRoboticsBook
using BenchmarkTools

robot = load_robot(:ur5)

q = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
v = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
a = [1.0, 1.5, 2.0, 0.5, 0.3, 0.1]
tau = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

println("=== ModernRoboticsBook.jl — UR5 (6-DOF) ===\n")

print("FK (space):              ")
display(@benchmark forward_kinematics_space($robot, $q))
println()

print("Jacobian (space):        ")
display(@benchmark jacobian_space($robot, $q))
println()

print("Inverse dynamics (RNEA): ")
display(@benchmark inverse_dynamics_rnea($robot, $q, $v, $a))
println()

print("Mass matrix (CRBA):      ")
display(@benchmark mass_matrix_crba($robot, $q))
println()

print("Mass matrix (n×RNEA):    ")
display(
    @benchmark mass_matrix_rnea(
        $q,
        $(robot.link_frames),
        $(robot.spatial_inertias),
        $(robot.screw_axes_space),
    )
)
println()

print("Gravity forces:          ")
display(@benchmark gravity_forces($robot, $q))
println()

print("Forward dynamics (CRBA): ")
display(@benchmark forward_dynamics_crba($robot, $q, $v, $tau))
println()

print("Forward dynamics (RNEA): ")
display(
    @benchmark forward_dynamics_rnea(
        $q,
        $v,
        $tau,
        $(robot.gravity),
        $(zeros(6)),
        $(robot.link_frames),
        $(robot.spatial_inertias),
        $(robot.screw_axes_space),
    )
)
println()

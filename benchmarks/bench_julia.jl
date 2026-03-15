#!/usr/bin/env julia
"""Benchmark ModernRoboticsBook.jl on UR5 for comparison with Pinocchio.

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

println("=== ModernRoboticsBook.jl — UR5 (6-DOF) ===\n")

print("FK (space):          ")
display(@benchmark forward_kinematics_space($robot, $q))
println()

print("FK (body):           ")
display(@benchmark forward_kinematics_body($robot, $q))
println()

print("Jacobian (space):    ")
display(@benchmark jacobian_space($robot, $q))
println()

print("Jacobian (body):     ")
display(@benchmark jacobian_body($robot, $q))
println()

print("Inverse dynamics:    ")
display(@benchmark inverse_dynamics($robot, $q, $v, $a))
println()

print("Mass matrix:         ")
display(@benchmark mass_matrix($robot, $q))
println()

print("Gravity forces:      ")
display(@benchmark gravity_forces($robot, $q))
println()

print("Forward dynamics:    ")
display(@benchmark forward_dynamics($robot, $q, $v, $a))
println()

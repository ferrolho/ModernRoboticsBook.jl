#!/usr/bin/env julia
"""Benchmark RigidBodyDynamics.jl on UR5 for comparison with ModernRoboticsBook.jl.

Run from the repo root:
    julia benchmarks/bench_rbd.jl

Requires RigidBodyDynamics.jl and BenchmarkTools.jl:
    julia -e 'using Pkg; Pkg.add(["RigidBodyDynamics", "BenchmarkTools"])'

The UR5 URDF is fetched via the robot_descriptions Python package:
    pip install robot-descriptions
    python -c 'from robot_descriptions import ur5_description'
"""

using RigidBodyDynamics
using BenchmarkTools

urdf_path = joinpath(
    homedir(),
    ".cache",
    "robot_descriptions",
    "example-robot-data",
    "robots",
    "ur_description",
    "urdf",
    "ur5_robot.urdf",
)

if !isfile(urdf_path)
    error(
        "UR5 URDF not found at $urdf_path. " *
        "Install it via: pip install robot-descriptions && " *
        "python -c 'from robot_descriptions import ur5_description'",
    )
end

mechanism = parse_urdf(urdf_path)
state = MechanismState(mechanism)

q = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
v = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]

set_configuration!(state, q)
set_velocity!(state, v)

# Build a SegmentedVector for accelerations (required by RBD's inverse_dynamics)
vd = similar(velocity(state))
vd .= [1.0, 1.5, 2.0, 0.5, 0.3, 0.1]

println("=== RigidBodyDynamics.jl — UR5 (6-DOF) ===\n")

print("Mass matrix:         ")
display(@benchmark mass_matrix($state))
println()

print("Mass matrix (0-alloc): ")
M = mass_matrix(state)
display(@benchmark mass_matrix!($M, $state))
println()

print("Inverse dynamics:    ")
display(@benchmark inverse_dynamics($state, $vd))
println()

print("Dynamics bias (c+g): ")
display(@benchmark dynamics_bias($state))
println()

print("Dynamics bias (0-alloc): ")
result = DynamicsResult(mechanism)
display(@benchmark dynamics_bias!($result, $state))
println()

# Performance Benchmarks

This page compares the performance of ModernRoboticsBook.jl against [Pinocchio](https://github.com/stack-of-tasks/pinocchio), a highly optimized C++ rigid-body dynamics library with Python bindings. All benchmarks use the **UR5** (6-DOF) robot at the same configuration.

## Results

| Function | ModernRoboticsBook.jl | Pinocchio 3.9 (C++) | Ratio |
|----------|----------------------:|--------------------:|------:|
| Forward kinematics | 6.7 μs | 0.86 μs | ~8x |
| Jacobian | 8.3 μs | 1.0 μs | ~8x |
| Inverse dynamics | 28.8 μs | 0.64 μs | ~45x |
| Mass matrix | 171 μs | 0.65 μs | ~263x |
| Gravity forces | 28.7 μs | 0.48 μs | ~60x |
| Forward dynamics | 258 μs | 2.82 μs | ~91x |

*Measured on Apple M2 (16 GB), Julia 1.12, Python 3.13. Julia timings are median values from BenchmarkTools.jl.*

## Analysis

The performance gap is expected. ModernRoboticsBook.jl is an **educational library** that follows the textbook algorithms directly for clarity, while Pinocchio is a **production library** written in templated C++ with hand-optimized memory layout.

The main factors behind the difference:

### Allocations

Our functions allocate temporary matrices on every call. For example, forward kinematics makes **524 allocations** per call, while Pinocchio pre-allocates all workspace in a `Data` struct and performs **zero allocations** per call.

| Function | Allocations | Memory |
|----------|------------:|-------:|
| Forward kinematics | 524 | 31 KiB |
| Jacobian | 647 | 39 KiB |
| Inverse dynamics | 2,304 | 135 KiB |
| Mass matrix | 13,862 | 812 KiB |
| Forward dynamics | 20,804 | 1.2 MiB |

### Algorithms

- **Mass matrix**: ModernRoboticsBook.jl calls `inverse_dynamics` *n* times with unit accelerations (the textbook approach). Pinocchio uses the **Composite Rigid Body Algorithm (CRBA)**, which computes the full matrix in a single pass.
- **Forward dynamics**: ModernRoboticsBook.jl explicitly forms and inverts the mass matrix. Pinocchio uses the **Articulated Body Algorithm (ABA)**, which solves for joint accelerations in O(n) without forming the mass matrix.

### When does this matter?

For **learning and prototyping**, ModernRoboticsBook.jl is fast enough — a full forward dynamics call takes ~258 μs, allowing ~3,800 evaluations per second. This is sufficient for trajectory optimization, offline simulation, and interactive exploration.

For **real-time control loops** (1 kHz+) or **large-scale optimization** (millions of evaluations), use a production library like Pinocchio, [RigidBodyDynamics.jl](https://github.com/JuliaRobotics/RigidBodyDynamics.jl), or [MuJoCo](https://mujoco.org/).

## Reproducing

```bash
# Julia benchmark
julia --project benchmarks/bench_julia.jl

# Pinocchio benchmark (requires: pip install pin robot-descriptions)
python benchmarks/bench_pinocchio.py
```

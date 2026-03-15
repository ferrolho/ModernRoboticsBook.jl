# Performance Benchmarks

This page compares the performance of ModernRoboticsBook.jl against [Pinocchio](https://github.com/stack-of-tasks/pinocchio), a highly optimized C++ rigid-body dynamics library with Python bindings. All benchmarks use the **UR5** (6-DOF) robot at the same configuration.

## Results

| Function | ModernRoboticsBook.jl | Pinocchio 3.9 (C++) | Ratio |
|----------|----------------------:|--------------------:|------:|
| Forward kinematics | 0.40 μs | 0.86 μs | **0.5x (faster)** |
| Jacobian | 0.45 μs | 1.0 μs | **0.5x (faster)** |
| Inverse dynamics | 3.64 μs | 0.64 μs | ~6x |
| Mass matrix | 22.0 μs | 0.65 μs | ~34x |
| Gravity forces | 3.68 μs | 0.48 μs | ~8x |
| Forward dynamics | 34.0 μs | 2.82 μs | ~12x |

*Measured on Apple M2 (16 GB), Julia 1.12, Python 3.13. Julia timings are median values from BenchmarkTools.jl.*

## Analysis

Forward kinematics and Jacobian computation are **faster than Pinocchio** thanks to Julia's StaticArrays, which eliminate heap allocations for small fixed-size matrices (4×4, 6×6).

The remaining gap in dynamics functions comes from algorithmic differences:

### Algorithms

- **Mass matrix**: ModernRoboticsBook.jl calls `inverse_dynamics` *n* times with unit accelerations (the textbook approach). Pinocchio uses the **Composite Rigid Body Algorithm (CRBA)**, which computes the full matrix in a single pass.
- **Forward dynamics**: ModernRoboticsBook.jl explicitly forms and inverts the mass matrix. Pinocchio uses the **Articulated Body Algorithm (ABA)**, which solves for joint accelerations in O(n) without forming the mass matrix.

### Allocations

With StaticArrays, allocations are dramatically reduced compared to the initial implementation:

| Function | Allocations | Memory |
|----------|------------:|-------:|
| Forward kinematics | 2 | 208 B |
| Jacobian | 2 | 384 B |
| Inverse dynamics | 186 | 13.4 KiB |
| Mass matrix | 1,114 | 80.7 KiB |
| Forward dynamics | 1,702 | 125.6 KiB |

### When does this matter?

For **learning and prototyping**, ModernRoboticsBook.jl is fast enough — a full forward dynamics call takes ~34 μs, allowing ~29,000 evaluations per second. This is sufficient for trajectory optimization, offline simulation, and interactive exploration.

For **real-time control loops** (1 kHz+) or **large-scale optimization** (millions of evaluations), use a production library like Pinocchio, [RigidBodyDynamics.jl](https://github.com/JuliaRobotics/RigidBodyDynamics.jl), or [MuJoCo](https://mujoco.org/).

## Reproducing

```bash
# Julia benchmark
julia --project benchmarks/bench_julia.jl

# Pinocchio benchmark (requires: pip install pin robot-descriptions)
python benchmarks/bench_pinocchio.py
```

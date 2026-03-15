# Performance Benchmarks

This page compares the performance of ModernRoboticsBook.jl against [Pinocchio](https://github.com/stack-of-tasks/pinocchio), a highly optimized C++ rigid-body dynamics library with Python bindings. All benchmarks use the **UR5** (6-DOF) robot at the same configuration.

## Results

| Function | ModernRoboticsBook.jl | In-place (`!`) | Pinocchio 3.9 (C++) | vs Pinocchio |
|----------|----------------------:|---------------:|--------------------:|-------------:|
| Forward kinematics | 0.31 μs | **0.30 μs** | 0.86 μs | **~3x faster** |
| Jacobian | 0.39 μs | **0.36 μs** | 1.0 μs | **~3x faster** |
| Inverse dynamics (RNEA) | 1.15 μs | **1.0 μs** | 0.64 μs | ~1.6x slower |
| Mass matrix (CRBA) | 1.14 μs | **0.91 μs** | 0.65 μs | ~1.4x slower |
| Gravity forces | 1.19 μs | — | 0.48 μs | ~2.5x slower |
| Forward dynamics (CRBA + RNEA) | 2.71 μs | — | 2.82 μs | **~1x (on par)** |

For reference, here are the timings of the textbook algorithms that were replaced by the optimized versions above:

| Function (textbook algorithm) | ModernRoboticsBook.jl | Pinocchio 3.9 (C++) | vs Pinocchio |
|-------------------------------|----------------------:|--------------------:|-------------:|
| Inverse dynamics (RNEA, unoptimized) | 3.81 μs | 0.64 μs | ~6x slower |
| Mass matrix (n × RNEA) | 22.2 μs | 0.65 μs | ~34x slower |
| Forward dynamics (M⁻¹ × ...) | 34.2 μs | 2.82 μs | ~12x slower |

*Measured on Apple M2 (16 GB), Julia 1.12, Python 3.13. Julia timings are median values from BenchmarkTools.jl.*

## Analysis

Forward kinematics and Jacobian computation are **~3x faster than Pinocchio** thanks to Julia's StaticArrays, which eliminate heap allocations for small fixed-size matrices (4×4, 6×6). The in-place variants (`forward_kinematics_space!`, `jacobian_space!`, etc.) achieve **zero allocations**.

Dynamics functions are now **within 1.4–2.5x of Pinocchio**, with forward dynamics on par. Both libraries use the same algorithms (RNEA, CRBA); the remaining gap is due to Pinocchio's hand-tuned C++/Eigen spatial algebra vs our general Julia `SMatrix` operations.

### Algorithms

- **Inverse dynamics**: Uses the **Recursive Newton-Euler Algorithm (RNEA)**, the same O(n) algorithm as Pinocchio. The in-place variant achieves zero allocations.
- **Mass matrix**: Uses the **Composite Rigid Body Algorithm (CRBA)**, which computes the full matrix in a single backward pass over composite spatial inertias — the same algorithm as Pinocchio. The in-place variant achieves zero allocations.
- **Forward dynamics**: Computes `M \ (τ - RNEA(q, dq, 0, g, F))` using CRBA for the mass matrix and a single RNEA call for the bias forces. Pinocchio uses the **Articulated Body Algorithm (ABA)**, which solves for joint accelerations in O(n) without forming the mass matrix.

### Allocations

| Function | Allocations | Memory | In-place |
|----------|------------:|-------:|---------:|
| Forward kinematics | 2 | 208 B | **0 (0 B)** |
| Jacobian | 2 | 384 B | **0 (0 B)** |
| Inverse dynamics (RNEA) | 13 | 3.4 KiB | **0 (0 B)** |
| Mass matrix (CRBA) | 9 | 4.6 KiB | **0 (0 B)** |
| Forward dynamics | 32 | 8.8 KiB | — |

### When does this matter?

For **learning and prototyping**, ModernRoboticsBook.jl is fast enough — a full forward dynamics call takes ~2.7 μs, allowing ~370,000 evaluations per second. This is sufficient for trajectory optimization, offline simulation, and interactive exploration.

For **real-time control loops** (1 kHz+) or **large-scale optimization** (millions of evaluations), the in-place variants bring performance close to production libraries like Pinocchio, [RigidBodyDynamics.jl](https://github.com/JuliaRobotics/RigidBodyDynamics.jl), or [MuJoCo](https://mujoco.org/).

## Reproducing

```bash
# Julia benchmark
julia --project benchmarks/bench_julia.jl

# Pinocchio benchmark (requires: pip install pin robot-descriptions)
python benchmarks/bench_pinocchio.py
```

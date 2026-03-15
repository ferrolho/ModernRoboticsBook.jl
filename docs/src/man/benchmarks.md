# Performance Benchmarks

This page compares the performance of ModernRoboticsBook.jl against [Pinocchio](https://github.com/stack-of-tasks/pinocchio) (C++) and [RigidBodyDynamics.jl](https://github.com/JuliaRobotics/RigidBodyDynamics.jl) (Julia). All benchmarks use the **UR5** (6-DOF) robot at the same configuration.

## Results

| Function | MRB.jl | MRB.jl (in-place) | Pinocchio 3.9 (C++) | RBD.jl |
|----------|-------:|---------:|--------------------:|-------:|
| Forward kinematics | 0.31 μs | **0.30 μs** | 0.86 μs | 0.41 μs |
| Jacobian | 0.39 μs | **0.36 μs** | 1.0 μs | 0.13 μs |
| Inverse dynamics (RNEA) | 1.15 μs | **1.0 μs** | 0.64 μs | 0.62 μs |
| Mass matrix (CRBA) | 1.14 μs | **0.91 μs** | 0.65 μs | 0.14 μs |
| Mass matrix (RNEA) | 6.79 μs | — | — | — |
| Gravity / dynamics bias | 1.19 μs | — | 0.48 μs | 0.55 μs |
| Forward dynamics (ABA) | 1.65 μs | **1.41 μs** | 2.82 μs | — |
| Forward dynamics (CRBA) | 2.80 μs | — | — | — |
| Forward dynamics (RNEA) | 5.44 μs | — | — | — |

*Measured on Apple M2 (16 GB), Julia 1.12, Python 3.13. Julia timings are median values from BenchmarkTools.jl. RBD.jl: the ee\_link frame is grabbed before removing fixed joints; timings are for in-place variants where available.*

## Analysis

Forward kinematics and Jacobian computation are **~3x faster than Pinocchio** thanks to Julia's StaticArrays, which eliminate heap allocations for small fixed-size matrices (4×4, 6×6). The in-place variants (`forward_kinematics_space!`, `jacobian_space!`, etc.) achieve **zero allocations**.

Dynamics functions are now **within 1.4–2x of Pinocchio**, with forward dynamics on par. Both libraries use the same algorithms (RNEA, CRBA); the remaining gap is due to Pinocchio's hand-tuned C++/Eigen spatial algebra vs our general Julia `SMatrix` operations.

RigidBodyDynamics.jl achieves extremely fast mass matrix (0.14 μs) and Jacobian (0.13 μs) computation through aggressive caching and code generation specialized to the mechanism topology. It uses a state-based API where `MechanismState` caches intermediate kinematic results, amortizing costs across queries. With fixed joints removed (the default), the kinematic tree is smaller and dynamics are faster. The `ee_link` frame is preserved on the parent body for FK queries even after pruning.

The textbook algorithm variants (`mass_matrix_rnea`, `forward_dynamics_rnea`) are provided for educational purposes — they directly mirror the textbook equations but are slower than the optimized defaults.

### Algorithms

- **Inverse dynamics**: Uses the **Recursive Newton-Euler Algorithm (RNEA)**, the same O(n) algorithm as Pinocchio and RBD.jl. The in-place variant achieves zero allocations.
- **Mass matrix**: Uses the **Composite Rigid Body Algorithm (CRBA)**, which computes the full matrix in a single backward pass over composite spatial inertias — the same algorithm as Pinocchio and RBD.jl. The in-place variant achieves zero allocations. The textbook variant `mass_matrix_rnea` calls `inverse_dynamics_rnea` n times with unit accelerations.
- **Forward dynamics**: Three implementations are available: `forward_dynamics_aba` uses the **Articulated Body Algorithm** — the same O(n) algorithm as Pinocchio — which solves for joint accelerations without forming the mass matrix. `forward_dynamics_crba` uses CRBA + a single RNEA call + backslash solve. `forward_dynamics_rnea` is the textbook variant that explicitly forms M⁻¹ and calls RNEA separately for each term.

### Allocations

| Function | Allocations | Memory | In-place |
|----------|------------:|-------:|---------:|
| Forward kinematics | 2 | 208 B | **0 (0 B)** |
| Jacobian | 2 | 384 B | **0 (0 B)** |
| Inverse dynamics (RNEA) | 13 | 3.4 KiB | **0 (0 B)** |
| Mass matrix (CRBA) | 9 | 4.6 KiB | **0 (0 B)** |
| Mass matrix (RNEA) | 76 | 20.7 KiB | — |
| Gravity forces | 17 | 3.6 KiB | — |
| Forward dynamics (ABA) | 23 | 6.1 KiB | **0 (0 B)** |
| Forward dynamics (CRBA) | 32 | 8.8 KiB | — |
| Forward dynamics (RNEA) | 78 | 19.6 KiB | — |

### When does this matter?

For **learning and prototyping**, ModernRoboticsBook.jl is fast enough — a full forward dynamics call takes ~2.7 μs, allowing ~370,000 evaluations per second. This is sufficient for trajectory optimization, offline simulation, and interactive exploration.

For **real-time control loops** (1 kHz+) or **large-scale optimization** (millions of evaluations), the in-place variants bring performance close to production libraries like Pinocchio and RigidBodyDynamics.jl.

## Reproducing

```bash
# ModernRoboticsBook.jl benchmark
julia --project benchmarks/bench_julia.jl

# RigidBodyDynamics.jl benchmark
julia benchmarks/bench_rbd.jl

# Pinocchio benchmark (requires: pip install pin robot-descriptions)
python benchmarks/bench_pinocchio.py
```

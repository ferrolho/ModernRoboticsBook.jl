# ModernRoboticsBook.jl

Julia port of the [official Modern Robotics library](https://github.com/NxRLab/ModernRobotics), which provides Python, MATLAB, and Mathematica implementations of the algorithms from [*Modern Robotics: Mechanics, Planning, and Control*](http://hades.mech.northwestern.edu/index.php/Modern_Robotics) by Kevin Lynch and Frank Park.

## Scope and Design

This package implements the algorithms from the Modern Robotics textbook for **serial chain robots** using the Product of Exponentials (PoE) formulation. It provides:

- **Simple, array-based API** — functions take plain arrays and return plain arrays. No complex type system or compilation step required.
- **Educational and optimized variants side by side** — textbook algorithms (e.g., `mass_matrix_rnea`) are available alongside production algorithms (e.g., `mass_matrix_crba`) with educational notes explaining the differences.
- **Competitive performance** — faster than Pinocchio (C++) for FK, Jacobian, and forward dynamics (ABA), with zero-allocation in-place variants for all core functions. See [Benchmarks](man/benchmarks.md) for details.
- **Robot model loading** — load robot models from JSON files (converted from URDF) via [`load_robot`](@ref).

**What this package is not**: a full-featured rigid body dynamics engine. It supports serial chains only (no branching kinematic trees, no closed-loop constraints, no collision geometry). For humanoids, quadrupeds, or other branching mechanisms, consider [RigidBodyDynamics.jl](https://github.com/JuliaRobotics/RigidBodyDynamics.jl), [Pinocchio](https://github.com/stack-of-tasks/pinocchio), or [MuJoCo](https://mujoco.org/).

## Installation

```julia
julia> import Pkg; Pkg.add("ModernRoboticsBook")
```

## API Reference

```@contents
Pages = [
    "lib/rigid-body-motions.md",
    "lib/forward-kinematics.md",
    "lib/velocity-kinematics.md",
    "lib/inverse-kinematics.md",
    "lib/dynamics.md",
    "lib/trajectory-generation.md",
    "lib/robot-control.md",
]
Depth = 1
```

## [Index](@id main-index)

```@index
Pages = [
    "lib/rigid-body-motions.md",
    "lib/forward-kinematics.md",
    "lib/velocity-kinematics.md",
    "lib/inverse-kinematics.md",
    "lib/dynamics.md",
    "lib/trajectory-generation.md",
    "lib/robot-control.md",
]
```

# ModernRoboticsBook.jl

[![Stable][docs-stable-img]][docs-stable-url]
[![Dev][docs-dev-img]][docs-dev-url]
[![CI][ci-img]][ci-url]
[![Codecov][codecov-img]][codecov-url]
[![Aqua QA][aqua-img]][aqua-url]

Julia port of the [official Modern Robotics library](https://github.com/NxRLab/ModernRobotics), which provides Python, MATLAB, and Mathematica implementations of the algorithms from [*Modern Robotics: Mechanics, Planning, and Control*](http://hades.mech.northwestern.edu/index.php/Modern_Robotics) by Kevin Lynch and Frank Park.

The package provides functions for rigid-body motions, forward/inverse kinematics, velocity kinematics, dynamics, trajectory generation, and robot control.

## Installation

```julia
julia> import Pkg; Pkg.add("ModernRoboticsBook")
```

## Quick start

```julia
using ModernRoboticsBook
import LinearAlgebra as LA

# Forward kinematics (body frame)
home_config = [
    -1  0  0  0
     0  1  0  6
     0  0 -1  2
     0  0  0  1.0
]
body_screw_axes = [0 0 -1 2 0   0
                   0 0  0 0 1   0
                   0 0  1 0 0 0.1]'
joint_positions = [π/2, 3, π]

T = fkin_body(home_config, body_screw_axes, joint_positions)
```

## Documentation

- [**Stable docs**](https://ferrolho.github.io/ModernRoboticsBook.jl/stable) — latest released version
- [**Dev docs**](https://ferrolho.github.io/ModernRoboticsBook.jl/dev) — current master branch
- [Examples](https://ferrolho.github.io/ModernRoboticsBook.jl/dev/man/examples/)
- [Function index](https://ferrolho.github.io/ModernRoboticsBook.jl/dev/#main-index-1)

## Attribution

This package is a Julia port of the [Modern Robotics library](https://github.com/NxRLab/ModernRobotics) accompanying the textbook:

> Kevin M. Lynch and Frank C. Park, *Modern Robotics: Mechanics, Planning, and Control*, Cambridge University Press, 2017, ISBN 9781107156302.

The textbook and companion materials are available at http://modernrobotics.org.

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://ferrolho.github.io/ModernRoboticsBook.jl/stable

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://ferrolho.github.io/ModernRoboticsBook.jl/dev

[ci-img]: https://github.com/ferrolho/ModernRoboticsBook.jl/actions/workflows/CI.yml/badge.svg?branch=master
[ci-url]: https://github.com/ferrolho/ModernRoboticsBook.jl/actions/workflows/CI.yml

[codecov-img]: https://codecov.io/gh/ferrolho/ModernRoboticsBook.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/ferrolho/ModernRoboticsBook.jl

[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

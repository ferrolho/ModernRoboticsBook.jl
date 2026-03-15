# Robot Model

## Type and Loading

```@autodocs
Modules = [ModernRoboticsBook]
Pages = ["robot.jl"]
Filter = t -> nameof(t) in [:Robot, :load_robot]
```

## Convenience Wrappers

The following functions accept a [`Robot`](@ref) as their first argument, automatically supplying the robot's kinematic and dynamic parameters. See the corresponding chapter pages for algorithmic details.

```@autodocs
Modules = [ModernRoboticsBook]
Pages = ["robot.jl"]
Filter = t -> nameof(t) ∉ [:Robot, :load_robot]
```

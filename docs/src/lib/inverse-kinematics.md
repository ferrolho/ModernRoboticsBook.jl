# Inverse Kinematics

*Textbook Chapter 6*

!!! tip "Robot convenience wrappers"
    For simplified calls using a [`Robot`](@ref) model (e.g., `inverse_kinematics_body(robot, target, guess, ε_ω, ε_v)`), see [Robot Model](@ref).

```@autodocs
Modules = [ModernRoboticsBook]
Pages = ["ModernRoboticsBook.jl"]
Filter = t -> nameof(t) in [:inverse_kinematics_body, :inverse_kinematics_space]
```

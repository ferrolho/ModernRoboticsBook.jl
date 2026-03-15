# Velocity Kinematics and Statics

*Textbook Chapter 5*

!!! tip "Robot convenience wrappers"
    For simplified calls using a [`Robot`](@ref) model (e.g., `jacobian_space(robot, q)`), see [Robot Model](@ref).

```@autodocs
Modules = [ModernRoboticsBook]
Pages = ["ModernRoboticsBook.jl"]
Filter = t -> nameof(t) in [:jacobian_body, :jacobian_body!, :jacobian_space, :jacobian_space!]
```

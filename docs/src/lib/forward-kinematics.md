# Forward Kinematics

*Textbook Chapter 4*

!!! tip "Robot convenience wrappers"
    For simplified calls using a [`Robot`](@ref) model (e.g., `forward_kinematics_space(robot, q)`), see [Robot Model](@ref).

```@autodocs
Modules = [ModernRoboticsBook]
Pages = ["ModernRoboticsBook.jl"]
Filter = t -> nameof(t) in [:forward_kinematics_body, :forward_kinematics_body!, :forward_kinematics_space, :forward_kinematics_space!]
```

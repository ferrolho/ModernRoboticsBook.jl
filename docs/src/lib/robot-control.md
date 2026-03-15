# Robot Control

*Textbook Chapter 11*

!!! tip "Robot convenience wrappers"
    For simplified calls using a [`Robot`](@ref) model (e.g., `computed_torque(robot, q, dq, ...)`), see [Robot Model](@ref).

```@autodocs
Modules = [ModernRoboticsBook]
Pages = ["ModernRoboticsBook.jl"]
Filter = t -> nameof(t) in [:computed_torque, :simulate_control]
```

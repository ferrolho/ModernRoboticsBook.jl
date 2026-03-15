# Dynamics of Open Chains

*Textbook Chapter 8*

!!! tip "Robot convenience wrappers"
    For simplified calls using a [`Robot`](@ref) model (e.g., `inverse_dynamics(robot, q, dq, ddq)`), see [Robot Model](@ref).

```@autodocs
Modules = [ModernRoboticsBook]
Pages = ["ModernRoboticsBook.jl"]
Filter = t -> nameof(t) in [
    :ad, :inverse_dynamics, :inverse_dynamics!, :mass_matrix, :mass_matrix!,
    :velocity_quadratic_forces, :gravity_forces, :end_effector_forces,
    :forward_dynamics, :euler_step,
    :inverse_dynamics_trajectory, :forward_dynamics_trajectory,
]
```

# Dynamics of Open Chains

*Textbook Chapter 8*

!!! tip "Robot convenience wrappers"
    For simplified calls using a [`Robot`](@ref) model (e.g., `inverse_dynamics_rnea(robot, q, dq, ddq)`), see [Robot Model](@ref).

```@autodocs
Modules = [ModernRoboticsBook]
Pages = ["ModernRoboticsBook.jl"]
Filter = t -> nameof(t) in [
    :ad, :inverse_dynamics_rnea, :inverse_dynamics_rnea!, :mass_matrix_crba, :mass_matrix_crba!,
    :mass_matrix_rnea, :velocity_quadratic_forces, :gravity_forces,
    :end_effector_forces, :forward_dynamics_aba, :forward_dynamics_aba!,
    :forward_dynamics_crba, :forward_dynamics_crba!,
    :forward_dynamics_rnea, :forward_dynamics_rnea!,
    :euler_step, :inverse_dynamics_trajectory, :forward_dynamics_trajectory,
]
```

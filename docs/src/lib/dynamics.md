# Dynamics of Open Chains

*Textbook Chapter 8*

```@autodocs
Modules = [ModernRoboticsBook]
Pages = ["ModernRoboticsBook.jl"]
Filter = t -> nameof(t) in [
    :ad, :inverse_dynamics, :mass_matrix, :velocity_quadratic_forces,
    :gravity_forces, :end_effector_forces, :forward_dynamics,
    :euler_step, :inverse_dynamics_trajectory, :forward_dynamics_trajectory,
]
```

# Trajectory Generation

*Textbook Chapter 9*

```@autodocs
Modules = [ModernRoboticsBook]
Pages = ["ModernRoboticsBook.jl"]
Filter = t -> nameof(t) in [
    :cubic_time_scaling, :quintic_time_scaling,
    :joint_trajectory, :screw_trajectory, :cartesian_trajectory,
]
```

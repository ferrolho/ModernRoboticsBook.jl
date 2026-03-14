# Rigid-Body Motions

*Textbook Chapter 3*

```@autodocs
Modules = [ModernRoboticsBook]
Pages = ["ModernRoboticsBook.jl"]
Filter = t -> nameof(t) in [
    :near_zero, :vec_to_so3, :so3_to_vec, :axis_angle3,
    :matrix_exp3, :matrix_log3, :rotation_position_to_transform, :transform_to_rotation_position, :transform_inv,
    :vec_to_se3, :se3_to_vec, :adjoint_representation, :screw_to_axis, :axis_angle6,
    :matrix_exp6, :matrix_log6, :project_to_so3, :project_to_se3,
    :distance_to_so3, :distance_to_se3, :is_so3, :is_se3,
]
```

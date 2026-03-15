export forward_kinematics_body,
    forward_kinematics_body!, forward_kinematics_space, forward_kinematics_space!

"""
    forward_kinematics_body(home_ee_pose, body_screw_axes, joint_positions)

Computes forward kinematics in the body frame for an open chain robot.

!!! info "Product of Exponentials (body frame)"
    Forward kinematics computes the end-effector pose from joint positions using the
    Product of Exponentials (PoE) formula. Each joint's contribution is a matrix
    exponential ``e^{[B_i]\\theta_i}`` that "moves" the end-effector from its home pose ``M``:

    ```math
    T(\\theta) = M \\, e^{[B_1]\\theta_1} \\, e^{[B_2]\\theta_2} \\cdots e^{[B_n]\\theta_n}
    ```

    The body-frame formulation expresses the screw axes as seen from the end-effector's
    perspective. This is equivalent to [`forward_kinematics_space`](@ref) (which uses
    ``T = e^{[S_1]\\theta_1} \\cdots e^{[S_n]\\theta_n} \\, M``) but can be more
    convenient when the task is defined relative to the tool.

# Arguments
- `home_ee_pose`: the ``4 \\times 4`` end-effector pose ``M \\in`` SE(3) when all joint positions are zero (home configuration).
- `body_screw_axes`: the joint screw axes ``B_i`` expressed in the end-effector (body) frame at the home configuration, as a ``6 \\times n`` matrix.
- `joint_positions`: an ``n``-vector of joint positions ``\\theta``.

# Returns
The ``4 \\times 4`` end-effector transformation matrix ``T \\in`` SE(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> home_ee_pose = [ -1  0  0  0 ;
                         0  1  0  6 ;
                         0  0 -1  2 ;
                         0  0  0  1 ];

julia> body_screw_axes = [  0  0 -1  2  0  0   ;
                            0  0  0  0  1  0   ;
                            0  0  1  0  0  0.1 ]';

julia> joint_positions = [ π/2, 3, π ];

julia> forward_kinematics_body(home_ee_pose, body_screw_axes, joint_positions)
4×4 Matrix{Float64}:
 -1.14424e-17  1.0           0.0  -5.0
  1.0          1.14424e-17   0.0   4.0
  0.0          0.0          -1.0   1.68584
  0.0          0.0           0.0   1.0
```
"""
function forward_kinematics_body(
    home_ee_pose::AbstractMatrix,
    body_screw_axes::AbstractMatrix,
    joint_positions::AbstractVector,
)
    T = similar(home_ee_pose, Float64, 4, 4)
    forward_kinematics_body!(T, home_ee_pose, body_screw_axes, joint_positions)
end

"""
    forward_kinematics_body!(T, home_ee_pose, body_screw_axes, joint_positions)

In-place version of [`forward_kinematics_body`](@ref). Writes the result into `T`.
"""
function forward_kinematics_body!(
    T::AbstractMatrix,
    home_ee_pose::AbstractMatrix,
    body_screw_axes::AbstractMatrix,
    joint_positions::AbstractVector,
)
    Ts = SMatrix{4,4}(home_ee_pose)
    for i in eachindex(joint_positions)
        Bi = SVector{6}(@view body_screw_axes[:, i])
        Ts = Ts * matrix_exp6(vec_to_se3(Bi * joint_positions[i]))
    end
    T .= Ts
    T
end

"""
    forward_kinematics_space(home_ee_pose, screw_axes, joint_positions)

Computes forward kinematics in the space frame for an open chain robot.

!!! info "Product of Exponentials (space frame)"
    Forward kinematics computes the end-effector pose from joint positions using the
    Product of Exponentials (PoE) formula. Each joint's contribution is a matrix
    exponential ``e^{[S_i]\\theta_i}`` applied to the home pose ``M``:

    ```math
    T(\\theta) = e^{[S_1]\\theta_1} \\, e^{[S_2]\\theta_2} \\cdots e^{[S_n]\\theta_n} \\, M
    ```

    The space-frame formulation expresses screw axes relative to a fixed base frame.
    This is equivalent to [`forward_kinematics_body`](@ref) (which uses
    ``T = M \\, e^{[B_1]\\theta_1} \\cdots e^{[B_n]\\theta_n}``) — both compute the
    same end-effector pose; the choice depends on which frame your screw axes are
    defined in.

# Arguments
- `home_ee_pose`: the ``4 \\times 4`` end-effector pose ``M \\in`` SE(3) when all joint positions are zero (home configuration).
- `screw_axes`: the joint screw axes ``S_i`` expressed in the space (base) frame at the home configuration, as a ``6 \\times n`` matrix.
- `joint_positions`: an ``n``-vector of joint positions ``\\theta``.

# Returns
The ``4 \\times 4`` end-effector transformation matrix ``T \\in`` SE(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> home_ee_pose = [ -1  0  0  0 ;
                         0  1  0  6 ;
                         0  0 -1  2 ;
                         0  0  0  1 ];

julia> screw_axes = [  0  0  1  4  0  0   ;
                       0  0  0  0  1  0   ;
                       0  0 -1 -6  0 -0.1 ]';

julia> joint_positions = [ π/2, 3, π ];

julia> forward_kinematics_space(home_ee_pose, screw_axes, joint_positions)
4×4 Matrix{Float64}:
 -1.14424e-17  1.0           0.0  -5.0
  1.0          1.14424e-17   0.0   4.0
  0.0          0.0          -1.0   1.68584
  0.0          0.0           0.0   1.0
```
"""
function forward_kinematics_space(
    home_ee_pose::AbstractMatrix,
    screw_axes::AbstractMatrix,
    joint_positions::AbstractVector,
)
    T = similar(home_ee_pose, Float64, 4, 4)
    forward_kinematics_space!(T, home_ee_pose, screw_axes, joint_positions)
end

"""
    forward_kinematics_space!(T, home_ee_pose, screw_axes, joint_positions)

In-place version of [`forward_kinematics_space`](@ref). Writes the result into `T`.
"""
function forward_kinematics_space!(
    T::AbstractMatrix,
    home_ee_pose::AbstractMatrix,
    screw_axes::AbstractMatrix,
    joint_positions::AbstractVector,
)
    Ts = SMatrix{4,4}(home_ee_pose)
    for i in reverse(eachindex(joint_positions))
        Si = SVector{6}(@view screw_axes[:, i])
        Ts = matrix_exp6(vec_to_se3(Si * joint_positions[i])) * Ts
    end
    T .= Ts
    T
end

export inverse_kinematics_body, inverse_kinematics_space

"""
    inverse_kinematics_body(body_screw_axes, home_ee_pose, target_config, initial_guess, angular_tolerance, linear_tolerance)

Computes inverse kinematics in the body frame for an open chain robot using Newton-Raphson iteration.

!!! info "Why is inverse kinematics iterative?"
    Unlike forward kinematics (which has a closed-form solution), inverse kinematics —
    finding joint angles for a desired end-effector pose — generally has no closed-form
    solution for arbitrary robots. This function uses Newton-Raphson iteration, repeatedly
    computing the error and updating the joint angles via the Jacobian until convergence.

# Arguments
- `body_screw_axes`: the joint screw axes ``B_i`` expressed in the end-effector (body) frame at the home configuration, as a ``6 \\times n`` matrix.
- `home_ee_pose`: the ``4 \\times 4`` end-effector pose ``M \\in`` SE(3) when all joint positions are zero (home configuration).
- `target_config`: the desired ``4 \\times 4`` end-effector pose ``T \\in`` SE(3).
- `initial_guess`: an ``n``-vector initial guess of joint positions ``\\theta_0``.
- `angular_tolerance`: small positive tolerance on the end-effector orientation error.
- `linear_tolerance`: small positive tolerance on the end-effector position error.

# Returns
A tuple `(joint_positions, success)` where `joint_positions` is the ``n``-vector solution and `success` is a `Bool`.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> body_screw_axes = [  0  0 -1  2  0  0   ;
                            0  0  0  0  1  0   ;
                            0  0  1  0  0  0.1 ]';

julia> home_ee_pose = [ -1  0  0  0 ;
                         0  1  0  6 ;
                         0  0 -1  2 ;
                         0  0  0  1 ];

julia> target_config = [ 0  1  0     -5 ;
                         1  0  0      4 ;
                         0  0 -1 1.6858 ;
                         0  0  0      1 ];

julia> initial_guess = [1.5, 2.5, 3];

julia> angular_tolerance, linear_tolerance = 0.01, 0.001;

julia> inverse_kinematics_body(body_screw_axes, home_ee_pose, target_config, initial_guess, angular_tolerance, linear_tolerance)
([1.5707381937148923, 2.999666997382943, 3.141539129217613], true)
```
"""
function inverse_kinematics_body(
    body_screw_axes::AbstractMatrix,
    home_ee_pose::AbstractMatrix,
    target_config::AbstractMatrix,
    initial_guess::AbstractVector,
    angular_tolerance::Number,
    linear_tolerance::Number,
)
    joint_positions = copy(initial_guess)
    i = 0
    maxiterations = 20
    Vb = se3_to_vec(
        matrix_log6(
            transform_inv(
                forward_kinematics_body(home_ee_pose, body_screw_axes, joint_positions),
            ) * target_config,
        ),
    )
    err = LA.norm(Vb[1:3]) > angular_tolerance || LA.norm(Vb[4:6]) > linear_tolerance
    while err && i < maxiterations
        joint_positions += LA.pinv(jacobian_body(body_screw_axes, joint_positions)) * Vb
        i += 1
        Vb = se3_to_vec(
            matrix_log6(
                transform_inv(
                    forward_kinematics_body(home_ee_pose, body_screw_axes, joint_positions),
                ) * target_config,
            ),
        )
        err = LA.norm(Vb[1:3]) > angular_tolerance || LA.norm(Vb[4:6]) > linear_tolerance
    end
    return joint_positions, !err
end

"""
    inverse_kinematics_space(screw_axes, home_ee_pose, target_config, initial_guess, angular_tolerance, linear_tolerance)

Computes inverse kinematics in the space frame for an open chain robot using Newton-Raphson iteration.

!!! info "Body vs space IK"
    This is the space-frame version of [`inverse_kinematics_body`](@ref). Both converge
    to the same solution; the difference is whether the error twist is computed in the
    body or space frame. The initial guess matters — the algorithm may converge to
    different solutions or fail to converge depending on the starting point.

# Arguments
- `screw_axes`: the joint screw axes ``S_i`` expressed in the space (base) frame at the home configuration, as a ``6 \\times n`` matrix.
- `home_ee_pose`: the ``4 \\times 4`` end-effector pose ``M \\in`` SE(3) when all joint positions are zero (home configuration).
- `target_config`: the desired ``4 \\times 4`` end-effector pose ``T \\in`` SE(3).
- `initial_guess`: an ``n``-vector initial guess of joint positions ``\\theta_0``.
- `angular_tolerance`: small positive tolerance on the end-effector orientation error.
- `linear_tolerance`: small positive tolerance on the end-effector position error.

# Returns
A tuple `(joint_positions, success)` where `joint_positions` is the ``n``-vector solution and `success` is a `Bool`.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> screw_axes = [  0  0  1  4  0  0   ;
                       0  0  0  0  1  0   ;
                       0  0 -1 -6  0 -0.1 ]';

julia> home_ee_pose = [ -1  0  0  0 ;
                         0  1  0  6 ;
                         0  0 -1  2 ;
                         0  0  0  1 ];

julia> target_config = [ 0  1  0     -5 ;
                         1  0  0      4 ;
                         0  0 -1 1.6858 ;
                         0  0  0      1 ];

julia> initial_guess = [1.5, 2.5, 3];

julia> angular_tolerance, linear_tolerance = 0.01, 0.001;

julia> inverse_kinematics_space(screw_axes, home_ee_pose, target_config, initial_guess, angular_tolerance, linear_tolerance)
([1.57073782965672, 2.999663844672525, 3.141534199856583], true)
```
"""
function inverse_kinematics_space(
    screw_axes::AbstractMatrix,
    home_ee_pose::AbstractMatrix,
    target_config::AbstractMatrix,
    initial_guess::AbstractVector,
    angular_tolerance::Number,
    linear_tolerance::Number,
)
    joint_positions = copy(initial_guess)
    i = 0
    maxiterations = 20
    Tsb = forward_kinematics_space(home_ee_pose, screw_axes, joint_positions)
    Vs =
        adjoint_representation(Tsb) *
        se3_to_vec(matrix_log6(transform_inv(Tsb) * target_config))
    err = LA.norm(Vs[1:3]) > angular_tolerance || LA.norm(Vs[4:6]) > linear_tolerance
    while err && i < maxiterations
        joint_positions += LA.pinv(jacobian_space(screw_axes, joint_positions)) * Vs
        i += 1
        Tsb = forward_kinematics_space(home_ee_pose, screw_axes, joint_positions)
        Vs =
            adjoint_representation(Tsb) *
            se3_to_vec(matrix_log6(transform_inv(Tsb) * target_config))
        err = LA.norm(Vs[1:3]) > angular_tolerance || LA.norm(Vs[4:6]) > linear_tolerance
    end
    joint_positions, !err
end

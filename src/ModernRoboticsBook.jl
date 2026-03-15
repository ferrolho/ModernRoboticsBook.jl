module ModernRoboticsBook

import LinearAlgebra as LA

export near_zero,
    vec_to_so3,
    so3_to_vec,
    axis_angle3,
    matrix_exp3,
    matrix_log3,
    rotation_position_to_transform,
    transform_to_rotation_position,
    transform_inv,
    vec_to_se3,
    se3_to_vec,
    adjoint_representation,
    screw_to_axis,
    axis_angle6,
    matrix_exp6,
    matrix_log6,
    project_to_so3,
    project_to_se3,
    distance_to_so3,
    distance_to_se3,
    is_so3,
    is_se3,
    forward_kinematics_body,
    forward_kinematics_space,
    jacobian_body,
    jacobian_space,
    inverse_kinematics_body,
    inverse_kinematics_space,
    ad,
    inverse_dynamics,
    mass_matrix,
    velocity_quadratic_forces,
    gravity_forces,
    end_effector_forces,
    forward_dynamics,
    euler_step,
    inverse_dynamics_trajectory,
    forward_dynamics_trajectory,
    cubic_time_scaling,
    quintic_time_scaling,
    joint_trajectory,
    screw_trajectory,
    cartesian_trajectory,
    computed_torque,
    simulate_control

# """
# *** BASIC HELPER FUNCTIONS ***
# """

"""
    near_zero(z)

Determines whether a scalar is small enough to be treated as zero.

# Arguments
- `z`: a scalar value.

# Returns
`true` if `z` is approximately zero (within absolute tolerance ``10^{-6}``); `false` otherwise.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> near_zero(-1e-7)
true
```
"""
near_zero(z::Number) = isapprox(z, zero(z); atol = 1e-6)

# """
# *** CHAPTER 3: RIGID-BODY MOTIONS ***
# """

"""
    vec_to_so3(ω)

Converts a 3-vector to an so(3) representation.

!!! info "What is so(3)?"
    so(3) is the space of all ``3 \\times 3`` skew-symmetric matrices. Each element
    represents an angular velocity. If ``\\omega`` is the angular velocity vector,
    then ``[\\omega]`` (its skew-symmetric form) can be used to compute cross products
    as matrix multiplications: ``\\omega \\times v = [\\omega] v``. This is the Lie
    algebra of the rotation group SO(3).

# Arguments
- `ω`: a 3-vector of angular velocities.

# Returns
The corresponding ``3 \\times 3`` skew-symmetric matrix in so(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> vec_to_so3([1, 2, 3])
3×3 Matrix{Int64}:
  0  -3   2
  3   0  -1
 -2   1   0
```
"""
function vec_to_so3(ω::AbstractVector)
    [
        0 -ω[3] ω[2]
        ω[3] 0 -ω[1]
        -ω[2] ω[1] 0
    ]
end

"""
    so3_to_vec(so3mat)

Converts an so(3) representation to a 3-vector.

# Arguments
- `so3mat`: a ``3 \\times 3`` skew-symmetric matrix in so(3).

# Returns
The corresponding 3-vector of angular velocities.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> so3_to_vec([0 -3 2; 3 0 -1; -2 1 0])
3-element Vector{Int64}:
 1
 2
 3
```
"""
function so3_to_vec(so3mat::AbstractMatrix)
    [so3mat[3, 2], so3mat[1, 3], so3mat[2, 1]]
end

"""
    axis_angle3(expc3)

Converts a 3-vector of exponential coordinates for rotation into axis-angle form.

# Arguments
- `expc3`: a 3-vector of exponential coordinates for rotation ``\\hat{\\omega}\\theta``.

# Returns
A tuple `(ω̂, θ)` of the unit rotation axis and the rotation angle.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> axis_angle3([1, 2, 3])
([0.2672612419124244, 0.5345224838248488, 0.8017837257372732], 3.7416573867739413)
```
"""
axis_angle3(expc3::AbstractVector) = LA.normalize(expc3), LA.norm(expc3)

"""
    matrix_exp3(so3mat)

Computes the matrix exponential of a matrix in so(3) using Rodrigues' formula.

!!! info "What is the matrix exponential?"
    The matrix exponential maps elements of a Lie algebra (infinitesimal motions)
    to the corresponding Lie group (finite motions). For rotations, it converts an
    angular velocity ``[\\omega]\\theta`` in so(3) into the actual rotation matrix
    ``R`` in SO(3) that results from rotating by angle ``\\theta`` about axis
    ``\\hat{\\omega}``. This is the mathematical foundation of the
    product-of-exponentials formula for robot kinematics.

# Arguments
- `so3mat`: a ``3 \\times 3`` skew-symmetric matrix in so(3).

# Returns
The corresponding rotation matrix ``R`` in SO(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> matrix_exp3([0 -3 2; 3 0 -1; -2 1 0])
3×3 Matrix{Float64}:
 -0.694921   0.713521  0.0892929
 -0.192007  -0.303785  0.933192 
  0.692978   0.63135   0.348107 
```
"""
function matrix_exp3(so3mat::AbstractMatrix)
    ωθ = so3_to_vec(so3mat)
    if near_zero(LA.norm(ωθ))
        return Matrix{Float64}(LA.I, 3, 3)
    else
        θ = axis_angle3(ωθ)[2]
        ωmat = so3mat / θ
        return LA.I + sin(θ) * ωmat + (1 - cos(θ)) * ωmat * ωmat
    end
end

"""
    matrix_log3(R)

Computes the matrix logarithm of a rotation matrix.

!!! info "What is the matrix logarithm?"
    The matrix logarithm is the inverse of the matrix exponential. It extracts
    the angular velocity ``[\\omega]\\theta`` from a rotation matrix ``R``,
    telling you the axis and angle of rotation that produced ``R``. This is
    essential for error computation in iterative algorithms like inverse kinematics.

# Arguments
- `R`: a rotation matrix in SO(3).

# Returns
The matrix logarithm of `R` in so(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> matrix_log3([0 0 1; 1 0 0; 0 1 0])
3×3 Matrix{Float64}:
  0.0     -1.2092   1.2092
  1.2092   0.0     -1.2092
 -1.2092   1.2092   0.0   
```
"""
function matrix_log3(R::AbstractMatrix)
    acosinput = (LA.tr(R) - 1) / 2
    if acosinput >= 1
        return zeros(3, 3)
    elseif acosinput <= -1
        if !near_zero(1 + R[3, 3])
            ω = (1 / √(2 * (1 + R[3, 3]))) * [R[1, 3], R[2, 3], 1 + R[3, 3]]
        elseif !near_zero(1 + R[2, 2])
            ω = (1 / √(2 * (1 + R[2, 2]))) * [R[1, 2], 1 + R[2, 2], R[3, 2]]
        else
            ω = (1 / √(2 * (1 + R[1, 1]))) * [1 + R[1, 1], R[2, 1], R[3, 1]]
        end
        return vec_to_so3(π * ω)
    else
        θ = acos(acosinput)
        return θ / 2 / sin(θ) * (R - R')
    end
end

"""
    rotation_position_to_transform(R, p)

Converts a rotation matrix and a position vector into homogeneous transformation matrix.

# Arguments
- `R`: a ``3 \\times 3`` rotation matrix.
- `p`: a 3-vector position.

# Returns
The corresponding ``4 \\times 4`` homogeneous transformation matrix ``T``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> rotation_position_to_transform([1 0 0; 0 0 -1; 0 1 0], [1, 2, 5])
4×4 Matrix{Int64}:
 1  0   0  1
 0  0  -1  2
 0  1   0  5
 0  0   0  1
```
"""
rotation_position_to_transform(R::AbstractMatrix, p::AbstractVector) =
    vcat(hcat(R, p), [0 0 0 1])

"""
    transform_to_rotation_position(T)

Converts a homogeneous transformation matrix into a rotation matrix and position vector.

# Arguments
- `T`: a ``4 \\times 4`` homogeneous transformation matrix.

# Returns
A tuple `(R, p)` of the rotation matrix and position vector.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> transform_to_rotation_position([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
([1 0 0; 0 0 -1; 0 1 0], [0, 0, 3])
```
"""
transform_to_rotation_position(T::AbstractMatrix) = T[1:3, 1:3], T[1:3, 4]

"""
    transform_inv(T)

Inverts a homogeneous transformation matrix.

# Arguments
- `T`: a ``4 \\times 4`` homogeneous transformation matrix.

# Returns
The inverse ``T^{-1}``, computed using ``R^T``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> transform_inv([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
4×4 Matrix{Int64}:
 1   0  0   0
 0   0  1  -3
 0  -1  0   0
 0   0  0   1
```
"""
function transform_inv(T::AbstractMatrix)
    R, p = transform_to_rotation_position(T)
    vcat(hcat(R', -R' * p), [0 0 0 1])
end

"""
    vec_to_se3(V)

Converts a spatial velocity vector into a 4x4 matrix in se(3).

!!! info "What is se(3)?"
    se(3) is the space of ``4 \\times 4`` matrices that represent twists — combined
    angular and linear velocities of a rigid body. It is the Lie algebra of the
    transformation group SE(3). A twist ``V = (\\omega, v)`` encodes both rotation
    (angular velocity ``\\omega``) and translation (linear velocity ``v``) in a
    single object, which is key to the product-of-exponentials formula for kinematics.

# Arguments
- `V`: a 6-vector spatial velocity (angular velocity, linear velocity).

# Returns
The corresponding ``4 \\times 4`` matrix in se(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> vec_to_se3([1, 2, 3, 4, 5, 6])
4×4 Matrix{Float64}:
  0.0  -3.0   2.0  4.0
  3.0   0.0  -1.0  5.0
 -2.0   1.0   0.0  6.0
  0.0   0.0   0.0  0.0
```
"""
vec_to_se3(V::AbstractVector) = vcat(hcat(vec_to_so3(V[1:3]), V[4:6]), zeros(1, 4))

"""
    se3_to_vec(se3mat)

Converts an se3 matrix into a spatial velocity vector.

# Arguments
- `se3mat`: a ``4 \\times 4`` matrix in se(3).

# Returns
The corresponding 6-vector twist (angular velocity, linear velocity).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> se3_to_vec([0 -3 2 4; 3 0 -1 5; -2 1 0 6; 0 0 0 0])
6-element Vector{Int64}:
 1
 2
 3
 4
 5
 6
```
"""
se3_to_vec(se3mat::AbstractMatrix) =
    vcat([se3mat[3, 2], se3mat[1, 3], se3mat[2, 1]], se3mat[1:3, 4])

"""
    adjoint_representation(T)

Computes the adjoint representation of a homogeneous transformation matrix.

!!! info "What is the adjoint representation?"
    The adjoint map transforms twists (spatial velocities) between reference frames.
    Given a twist ``V`` expressed in one frame, ``[\\text{Ad}_T] V`` gives the same
    physical motion expressed in the frame related by ``T``. This is used throughout
    robotics to transform Jacobians and wrenches between body and space frames.

# Arguments
- `T`: a ``4 \\times 4`` homogeneous transformation matrix.

# Returns
The ``6 \\times 6`` adjoint representation ``[\\text{Ad}_T]``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> adjoint_representation([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
6×6 Matrix{Float64}:
 1.0  0.0   0.0  0.0  0.0   0.0
 0.0  0.0  -1.0  0.0  0.0   0.0
 0.0  1.0   0.0  0.0  0.0   0.0
 0.0  0.0   3.0  1.0  0.0   0.0
 3.0  0.0   0.0  0.0  0.0  -1.0
 0.0  0.0   0.0  0.0  1.0   0.0
```
"""
function adjoint_representation(T::AbstractMatrix)
    R, p = transform_to_rotation_position(T)
    vcat(hcat(R, zeros(3, 3)), hcat(vec_to_so3(p) * R, R))
end

"""
    screw_to_axis(q, s, h)

Takes a parametric description of a screw axis and converts it to a normalized screw axis.

!!! info "What is a screw axis?"
    A screw axis describes any rigid-body motion as a rotation about an axis combined
    with a translation along that axis (like a corkscrew). It is defined by a point
    ``q`` on the axis, a direction ``s``, and a pitch ``h`` (ratio of linear to angular
    motion). When ``h = 0`` the motion is pure rotation; when ``h = \\infty`` it is pure
    translation. Every joint in a robot can be described by a screw axis.

# Arguments
- `q`: a point on the screw axis.
- `s`: the unit direction vector of the screw axis.
- `h`: the pitch of the screw axis.

# Returns
The normalized 6-vector screw axis ``\\mathcal{S}``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> screw_to_axis([3; 0; 0], [0; 0; 1], 2)
6-element Vector{Int64}:
  0
  0
  1
  0
 -3
  2
```
"""
screw_to_axis(q::AbstractVector, s::AbstractVector, h::Number) =
    vcat(s, LA.cross(q, s) + h * s)

"""
    axis_angle6(expc6)

Converts a 6-vector of exponential coordinates into screw axis-angle form.

# Arguments
- `expc6`: a 6-vector of exponential coordinates ``\\mathcal{S}\\theta``.

# Returns
A tuple `(S, θ)` of the screw axis and the distance travelled along the axis.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> axis_angle6([1, 0, 0, 1, 2, 3])
([1.0, 0.0, 0.0, 1.0, 2.0, 3.0], 1.0)
```
"""
function axis_angle6(expc6::AbstractVector)
    θ = LA.norm(expc6[1:3])
    if near_zero(θ)
        θ = LA.norm(expc6[4:6])
    end
    expc6 / θ, θ
end

"""
    matrix_exp6(se3mat)

Computes the matrix exponential of an se(3) representation of exponential coordinates.

!!! info "What does this do in SE(3)?"
    This is the SE(3) version of [`matrix_exp3`](@ref). It converts a twist
    ``[\\mathcal{S}]\\theta`` in se(3) into the rigid-body transformation ``T``
    in SE(3) that results from following that screw motion. This is how each
    joint's contribution is computed in the product-of-exponentials formula
    for forward kinematics.

# Arguments
- `se3mat`: a ``4 \\times 4`` matrix in se(3).

# Returns
The corresponding homogeneous transformation matrix ``T`` in SE(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> matrix_exp6([0 0 0 0; 0 0 -1.57079632 2.35619449; 0 1.57079632 0 2.35619449; 0 0 0 0])
4×4 Matrix{Float64}:
 1.0  0.0         0.0        0.0       
 0.0  6.7949e-9  -1.0        1.01923e-8
 0.0  1.0         6.7949e-9  3.0       
 0.0  0.0         0.0        1.0       
```
"""
function matrix_exp6(se3mat::AbstractMatrix)
    ωθ = so3_to_vec(se3mat[1:3, 1:3])
    if near_zero(LA.norm(ωθ))
        return vcat(hcat(Matrix{Float64}(LA.I, 3, 3), se3mat[1:3, 4]), [0 0 0 1])
    else
        θ = axis_angle3(ωθ)[2]
        ωmat = se3mat[1:3, 1:3] / θ
        return vcat(
            hcat(
                matrix_exp3(se3mat[1:3, 1:3]),
                (
                    Matrix{Float64}(LA.I, 3, 3) * θ +
                    (1 - cos(θ)) * ωmat +
                    (θ - sin(θ)) * ωmat * ωmat
                ) * se3mat[1:3, 4] / θ,
            ),
            [0 0 0 1],
        )
    end
end

"""
    matrix_log6(T)

Computes the matrix logarithm of a homogeneous transformation matrix.

!!! info "What does this do in SE(3)?"
    This is the SE(3) version of [`matrix_log3`](@ref). It extracts the twist
    ``[\\mathcal{S}]\\theta`` from a transformation ``T``, recovering the screw
    motion that produced it. This is used in inverse kinematics to compute the
    body twist error between the current and desired end-effector configurations.

# Arguments
- `T`: a ``4 \\times 4`` homogeneous transformation matrix in SE(3).

# Returns
The ``4 \\times 4`` matrix logarithm ``[\\mathcal{S}\\theta]`` in se(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> matrix_log6([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
4×4 Matrix{Float64}:
 0.0  0.0      0.0     0.0    
 0.0  0.0     -1.5708  2.35619
 0.0  1.5708   0.0     2.35619
 0.0  0.0      0.0     0.0    
```
"""
function matrix_log6(T::AbstractMatrix)
    R, p = transform_to_rotation_position(T)
    ωmat = matrix_log3(R)
    if iszero(ωmat)
        return vcat(hcat(zeros(3, 3), T[1:3, 4]), [0 0 0 0])
    else
        θ = acos((LA.tr(R) - 1) / 2)
        return vcat(
            hcat(
                ωmat,
                (LA.I - ωmat / 2 + (1 / θ - 1 / tan(θ / 2) / 2) * ωmat * ωmat / θ) *
                T[1:3, 4],
            ),
            [0 0 0 0],
        )
    end
end

"""
    project_to_so3(mat)

Returns a projection of mat into SO(3).

!!! info "Why project to SO(3)?"
    Numerical computations (e.g. repeated matrix multiplications) can cause rotation
    matrices to drift away from SO(3) — they may no longer be exactly orthogonal with
    determinant 1. This function finds the closest valid rotation matrix using SVD,
    which is essential for maintaining physical consistency in long-running simulations.

# Arguments
- `mat`: a ``3 \\times 3`` matrix near SO(3).

# Returns
The closest rotation matrix in SO(3), computed using SVD projection.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> project_to_so3([0.675 0.150  0.720; 0.370 0.771 -0.511; -0.630 0.619  0.472])
3×3 Matrix{Float64}:
  0.679011  0.148945   0.718859
  0.373207  0.773196  -0.512723
 -0.632187  0.616428   0.469421
```
"""
function project_to_so3(mat::AbstractMatrix)
    F = LA.svd(mat)
    R = F.U * F.Vt
    if LA.det(R) < 0
        R = F.U * LA.Diagonal([1, 1, -1]) * F.Vt
    end
    return R
end

"""
    project_to_se3(mat)

Returns a projection of mat into SE(3).

!!! info "Why project to SE(3)?"
    Similar to [`project_to_so3`](@ref), this corrects numerical drift in homogeneous
    transformation matrices. After many multiplications, the rotation part may no longer
    be orthogonal. This projects the matrix back onto SE(3) to maintain valid rigid-body
    transformations.

# Arguments
- `mat`: a ``4 \\times 4`` matrix near SE(3).

# Returns
The closest homogeneous transformation matrix in SE(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> project_to_se3([0.675 0.150 0.720 1.2; 0.370 0.771 -0.511 5.4; -0.630 0.619 0.472 3.6; 0.003 0.002 0.010 0.9])
4×4 Matrix{Float64}:
  0.679011  0.148945   0.718859  1.2
  0.373207  0.773196  -0.512723  5.4
 -0.632187  0.616428   0.469421  3.6
  0.0       0.0        0.0       1.0
```
"""
project_to_se3(mat::AbstractMatrix) =
    rotation_position_to_transform(project_to_so3(mat[1:3, 1:3]), mat[1:3, 4])

"""
    distance_to_so3(mat)

Returns the Frobenius norm to describe the distance of mat from the SO(3) manifold.

# Arguments
- `mat`: a ``3 \\times 3`` matrix.

# Returns
The Frobenius norm of ``R^T R - I``, or ``10^9`` if ``\\det(R) < 0``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> distance_to_so3([1.0 0.0 0.0; 0.0 0.1 -0.95; 0.0 1.0 0.1])
0.08835298523536149
```
"""
distance_to_so3(mat::AbstractMatrix) = LA.det(mat) > 0 ? LA.norm(mat'mat - LA.I) : 1e+9

"""
    distance_to_se3(mat)

Returns the Frobenius norm to describe the distance of mat from the SE(3) manifold.

# Arguments
- `mat`: a ``4 \\times 4`` matrix.

# Returns
The Frobenius distance from SE(3), or ``10^9`` if ``\\det(R) < 0``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> distance_to_se3([1.0 0.0 0.0 1.2; 0.0 0.1 -0.95 1.5; 0.0 1.0 0.1 -0.9; 0.0 0.0 0.1 0.98])
0.13493053768513638
```
"""
function distance_to_se3(mat::AbstractMatrix)
    matR = mat[1:3, 1:3]
    if LA.det(matR) > 0
        LA.norm(hcat(vcat(matR'matR, zeros(1, 3)), mat[4, :]) - LA.I)
    else
        1e+9
    end
end

"""
    is_so3(mat)

Returns true if mat is close to or on the manifold SO(3).

# Arguments
- `mat`: a ``3 \\times 3`` matrix.

# Returns
`true` if `mat` is close to SO(3); `false` otherwise.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> is_so3([1.0 0.0 0.0; 0.0 0.1 -0.95; 0.0 1.0 0.1])
false
```
"""
is_so3(mat::AbstractMatrix) = abs(distance_to_so3(mat)) < 1e-3

"""
    is_se3(mat)

Returns true if mat is close to or on the manifold SE(3).

# Arguments
- `mat`: a ``4 \\times 4`` matrix.

# Returns
`true` if `mat` is close to SE(3); `false` otherwise.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> is_se3([1.0 0.0 0.0 1.2; 0.0 0.1 -0.95 1.5; 0.0 1.0 0.1 -0.9; 0.0 0.0 0.1 0.98])
false
```
"""
is_se3(mat::AbstractMatrix) = abs(distance_to_se3(mat)) < 1e-3

# """
# *** CHAPTER 4: FORWARD KINEMATICS ***
# """

"""
    forward_kinematics_body(home_config, body_screw_axes, joint_positions)

Computes forward kinematics in the body frame for an open chain robot.

!!! info "What is the body frame?"
    In robotics, we often need to know where the end-effector is given the joint angles.
    The "body frame" formulation expresses the joint screw axes as seen from the
    end-effector's perspective. This is equivalent to `forward_kinematics_space` but can
    be more convenient when the task is defined relative to the tool.

# Arguments
- `home_config`: the ``4 \\times 4`` home configuration ``M`` of the end-effector (SE(3)).
- `body_screw_axes`: the joint screw axes ``B_i`` in the end-effector (body) frame, as a ``6 \\times n`` matrix.
- `joint_positions`: an ``n``-vector of joint positions ``\\theta``.

# Returns
The ``4 \\times 4`` end-effector transformation matrix ``T \\in`` SE(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> home_config = [ -1  0  0  0 ;
              0  1  0  6 ;
              0  0 -1  2 ;
              0  0  0  1 ];

julia> body_screw_axes = [  0  0 -1  2  0  0   ;
                  0  0  0  0  1  0   ;
                  0  0  1  0  0  0.1 ]';

julia> joint_positions = [ π/2, 3, π ];

julia> forward_kinematics_body(home_config, body_screw_axes, joint_positions)
4×4 Matrix{Float64}:
 -1.14424e-17  1.0           0.0  -5.0
  1.0          1.14424e-17   0.0   4.0
  0.0          0.0          -1.0   1.68584
  0.0          0.0           0.0   1.0
```
"""
function forward_kinematics_body(
    home_config::AbstractMatrix,
    body_screw_axes::AbstractMatrix,
    joint_positions::AbstractVector,
)
    for i in eachindex(joint_positions)
        home_config *= matrix_exp6(vec_to_se3(body_screw_axes[:, i] * joint_positions[i]))
    end
    home_config
end

"""
    forward_kinematics_space(home_config, screw_axes, joint_positions)

Computes forward kinematics in the space frame for an open chain robot.

!!! info "What is the space frame?"
    The "space frame" formulation expresses joint screw axes relative to a fixed base
    frame. Both body and space formulations compute the same end-effector pose — the
    choice depends on which frame your screw axes are defined in.

# Arguments
- `home_config`: the ``4 \\times 4`` home configuration ``M`` of the end-effector (SE(3)).
- `screw_axes`: the joint screw axes ``S_i`` in the space frame, as a ``6 \\times n`` matrix.
- `joint_positions`: an ``n``-vector of joint positions ``\\theta``.

# Returns
The ``4 \\times 4`` end-effector transformation matrix ``T \\in`` SE(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> home_config = [ -1  0  0  0 ;
              0  1  0  6 ;
              0  0 -1  2 ;
              0  0  0  1 ];

julia> screw_axes = [  0  0  1  4  0  0   ;
                  0  0  0  0  1  0   ;
                  0  0 -1 -6  0 -0.1 ]';

julia> joint_positions = [ π/2, 3, π ];

julia> forward_kinematics_space(home_config, screw_axes, joint_positions)
4×4 Matrix{Float64}:
 -1.14424e-17  1.0           0.0  -5.0
  1.0          1.14424e-17   0.0   4.0
  0.0          0.0          -1.0   1.68584
  0.0          0.0           0.0   1.0
```
"""
function forward_kinematics_space(
    home_config::AbstractMatrix,
    screw_axes::AbstractMatrix,
    joint_positions::AbstractVector,
)
    for i in reverse(eachindex(joint_positions))
        home_config =
            matrix_exp6(vec_to_se3(screw_axes[:, i] * joint_positions[i])) * home_config
    end
    home_config
end

# """
# *** CHAPTER 5: VELOCITY KINEMATICS AND STATICS ***
# """

"""
    jacobian_body(body_screw_axes, joint_positions)

Computes the body Jacobian ``J_b(\\theta)`` for an open chain robot.

!!! info "What is the Jacobian?"
    The Jacobian relates joint velocities to end-effector velocity:
    ``V_b = J_b(\\theta) \\dot{\\theta}``. It tells you how fast and in which direction
    the end-effector moves for a given set of joint velocities. It is also used to map
    wrenches (forces/torques) from the end-effector back to the joints.

# Arguments
- `body_screw_axes`: the joint screw axes ``B_i`` in the end-effector (body) frame at the home position, as a ``6 \\times n`` matrix.
- `joint_positions`: an ``n``-vector of joint positions ``\\theta``.

# Returns
The ``6 \\times n`` body Jacobian ``J_b(\\theta)``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> body_screw_axes = [0 0 1   0 0.2 0.2;
                1 0 0   2   0   3;
                0 1 0   0   2   1;
                1 0 0 0.2 0.3 0.4]';

julia> joint_positions = [0.2, 1.1, 0.1, 1.2];

julia> jacobian_body(body_screw_axes, joint_positions)
6×4 Matrix{Float64}:
 -0.0452841  0.995004    0.0       1.0
  0.743593   0.0930486   0.362358  0.0
 -0.667097   0.0361754  -0.932039  0.0
  2.32586    1.66809     0.564108  0.2
 -1.44321    2.94561     1.43307   0.3
 -2.0664     1.82882    -1.58869   0.4
```
"""
function jacobian_body(body_screw_axes::AbstractMatrix, joint_positions::AbstractVector)
    T = LA.I
    Jb = copy(body_screw_axes)
    for i in Iterators.reverse(firstindex(joint_positions):(lastindex(joint_positions)-1))
        T *= matrix_exp6(vec_to_se3(body_screw_axes[:, i+1] * -joint_positions[i+1]))
        Jb[:, i] = adjoint_representation(T) * body_screw_axes[:, i]
    end
    Jb
end

"""
    jacobian_space(screw_axes, joint_positions)

Computes the space Jacobian ``J_s(\\theta)`` for an open chain robot.

!!! info "Body vs space Jacobian"
    The space Jacobian ``J_s`` gives the end-effector velocity in the fixed space frame,
    while the body Jacobian ``J_b`` gives it in the end-effector's own frame. They are
    related by ``J_s = [\\text{Ad}_{T_{sb}}] J_b``. Use whichever matches the frame your
    task is defined in.

# Arguments
- `screw_axes`: the joint screw axes ``S_i`` in the space frame at the home position, as a ``6 \\times n`` matrix.
- `joint_positions`: an ``n``-vector of joint positions ``\\theta``.

# Returns
The ``6 \\times n`` space Jacobian ``J_s(\\theta)``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> screw_axes = [0 0 1   0 0.2 0.2;
                1 0 0   2   0   3;
                0 1 0   0   2   1;
                1 0 0 0.2 0.3 0.4]';

julia> joint_positions = [0.2, 1.1, 0.1, 1.2];

julia> jacobian_space(screw_axes, joint_positions)
6×4 Matrix{Float64}:
 0.0  0.980067  -0.0901156   0.957494 
 0.0  0.198669   0.444554    0.284876 
 1.0  0.0        0.891207   -0.0452841
 0.0  1.95219   -2.21635    -0.511615 
 0.2  0.436541  -2.43713     2.77536  
 0.2  2.96027    3.23573     2.22512  
```
"""
function jacobian_space(screw_axes::AbstractMatrix, joint_positions::AbstractVector)
    T = LA.I
    Js = copy(screw_axes)
    for i = (firstindex(joint_positions)+1):lastindex(joint_positions)
        T *= matrix_exp6(vec_to_se3(screw_axes[:, i-1] * joint_positions[i-1]))
        Js[:, i] = adjoint_representation(T) * screw_axes[:, i]
    end
    Js
end

# """
# *** CHAPTER 6: INVERSE KINEMATICS ***
# """

"""
    inverse_kinematics_body(body_screw_axes, home_config, target_config, initial_guess, angular_tolerance, linear_tolerance)

Computes inverse kinematics in the body frame for an open chain robot using Newton-Raphson iteration.

!!! info "Why is inverse kinematics iterative?"
    Unlike forward kinematics (which has a closed-form solution), inverse kinematics —
    finding joint angles for a desired end-effector pose — generally has no closed-form
    solution for arbitrary robots. This function uses Newton-Raphson iteration, repeatedly
    computing the error and updating the joint angles via the Jacobian until convergence.

# Arguments
- `body_screw_axes`: the joint screw axes ``B_i`` in the end-effector (body) frame, as a ``6 \\times n`` matrix.
- `home_config`: the ``4 \\times 4`` home configuration ``M`` of the end-effector (SE(3)).
- `target_config`: the desired ``4 \\times 4`` end-effector configuration ``T`` (SE(3)).
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

julia> home_config = [ -1  0  0  0 ;
              0  1  0  6 ;
              0  0 -1  2 ;
              0  0  0  1 ];

julia> target_config = [ 0  1  0     -5 ;
             1  0  0      4 ;
             0  0 -1 1.6858 ;
             0  0  0      1 ];

julia> initial_guess = [1.5, 2.5, 3];

julia> angular_tolerance, linear_tolerance = 0.01, 0.001;

julia> inverse_kinematics_body(body_screw_axes, home_config, target_config, initial_guess, angular_tolerance, linear_tolerance)
([1.5707381937148923, 2.999666997382943, 3.141539129217613], true)
```
"""
function inverse_kinematics_body(
    body_screw_axes::AbstractMatrix,
    home_config::AbstractMatrix,
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
                forward_kinematics_body(home_config, body_screw_axes, joint_positions),
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
                    forward_kinematics_body(home_config, body_screw_axes, joint_positions),
                ) * target_config,
            ),
        )
        err = LA.norm(Vb[1:3]) > angular_tolerance || LA.norm(Vb[4:6]) > linear_tolerance
    end
    return joint_positions, !err
end

"""
    inverse_kinematics_space(screw_axes, home_config, target_config, initial_guess, angular_tolerance, linear_tolerance)

Computes inverse kinematics in the space frame for an open chain robot using Newton-Raphson iteration.

!!! info "Body vs space IK"
    This is the space-frame version of [`inverse_kinematics_body`](@ref). Both converge
    to the same solution; the difference is whether the error twist is computed in the
    body or space frame. The initial guess matters — the algorithm may converge to
    different solutions or fail to converge depending on the starting point.

# Arguments
- `screw_axes`: the joint screw axes ``S_i`` in the space frame, as a ``6 \\times n`` matrix.
- `home_config`: the ``4 \\times 4`` home configuration ``M`` of the end-effector (SE(3)).
- `target_config`: the desired ``4 \\times 4`` end-effector configuration ``T`` (SE(3)).
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

julia> home_config = [ -1  0  0  0 ;
              0  1  0  6 ;
              0  0 -1  2 ;
              0  0  0  1 ];

julia> target_config = [ 0  1  0     -5 ;
             1  0  0      4 ;
             0  0 -1 1.6858 ;
             0  0  0      1 ];

julia> initial_guess = [1.5, 2.5, 3];

julia> angular_tolerance, linear_tolerance = 0.01, 0.001;

julia> inverse_kinematics_space(screw_axes, home_config, target_config, initial_guess, angular_tolerance, linear_tolerance)
([1.57073782965672, 2.999663844672525, 3.141534199856583], true)
```
"""
function inverse_kinematics_space(
    screw_axes::AbstractMatrix,
    home_config::AbstractMatrix,
    target_config::AbstractMatrix,
    initial_guess::AbstractVector,
    angular_tolerance::Number,
    linear_tolerance::Number,
)
    joint_positions = copy(initial_guess)
    i = 0
    maxiterations = 20
    Tsb = forward_kinematics_space(home_config, screw_axes, joint_positions)
    Vs =
        adjoint_representation(Tsb) *
        se3_to_vec(matrix_log6(transform_inv(Tsb) * target_config))
    err = LA.norm(Vs[1:3]) > angular_tolerance || LA.norm(Vs[4:6]) > linear_tolerance
    while err && i < maxiterations
        joint_positions += LA.pinv(jacobian_space(screw_axes, joint_positions)) * Vs
        i += 1
        Tsb = forward_kinematics_space(home_config, screw_axes, joint_positions)
        Vs =
            adjoint_representation(Tsb) *
            se3_to_vec(matrix_log6(transform_inv(Tsb) * target_config))
        err = LA.norm(Vs[1:3]) > angular_tolerance || LA.norm(Vs[4:6]) > linear_tolerance
    end
    joint_positions, !err
end

# """
# *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
# """

"""
    ad(V)

Computes the ``6 \\times 6`` matrix ``[\\text{ad}_V]`` used to calculate the Lie bracket ``[V_1, V_2] = [\\text{ad}_{V_1}] V_2``.

!!! info "What is the Lie bracket?"
    The Lie bracket measures how two twists interact — it captures the velocity
    produced by the non-commutativity of two motions. In dynamics, the ``[\\text{ad}]``
    matrix appears in the recursive Newton-Euler algorithm, where it accounts for
    how the velocity of one link affects the forces on adjacent links through
    Coriolis and centripetal effects.

# Arguments
- `V`: a 6-vector spatial velocity (twist).

# Returns
The ``6 \\times 6`` matrix ``[\\text{ad}_V]``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> ad([1, 2, 3, 4, 5, 6])
6×6 Matrix{Float64}:
  0.0  -3.0   2.0   0.0   0.0   0.0
  3.0   0.0  -1.0   0.0   0.0   0.0
 -2.0   1.0   0.0   0.0   0.0   0.0
  0.0  -6.0   5.0   0.0  -3.0   2.0
  6.0   0.0  -4.0   3.0   0.0  -1.0
 -5.0   4.0   0.0  -2.0   1.0   0.0
```
"""
function ad(V::AbstractVector)
    ωmat = vec_to_so3(V[1:3])
    vcat(hcat(ωmat, zeros(3, 3)), hcat(vec_to_so3(V[4:6]), ωmat))
end

"""
    inverse_dynamics(joint_positions, joint_velocities, joint_accelerations, gravity, tip_wrench, link_frames, spatial_inertias, screw_axes)

Computes inverse dynamics in the space frame for an open chain robot using
forward-backward Newton-Euler iterations.

!!! info "What is inverse dynamics?"
    Inverse dynamics answers: "what joint torques do I need to produce a desired motion?"
    Given the joint positions, velocities, and accelerations, it computes the required
    forces using the recursive Newton-Euler algorithm. This is the workhorse behind most
    dynamics computations — the mass matrix, Coriolis terms, and gravity vector are all
    computed by calling this function with special inputs.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `joint_velocities`: the ``n``-vector of joint rates.
- `joint_accelerations`: the ``n``-vector of joint accelerations.
- `gravity`: the gravity vector ``g`` (e.g., `[0, 0, -9.8]`).
- `tip_wrench`: the wrench ``\\mathcal{F}_{\\text{tip}}`` applied by the end-effector expressed in frame ``\\{n+1\\}``.
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}`` (the end-effector frame relative to the last link).
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.

# Returns
The ``n``-vector of required joint forces/torques.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> inverse_dynamics([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], [2, 1.5, 1], [0, 0, -9.8], [1, 1, 1, 1, 1, 1], link_frames, spatial_inertias, screw_axes)
3-element Vector{Float64}:
  74.6961615528745
 -33.067660158514585
  -3.2305731379014246
```
"""
function inverse_dynamics(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_accelerations::AbstractVector,
    gravity::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    Mi = LA.I
    Ai = zeros(eltype(joint_positions), 6, n)
    AdTi = Array{Array{eltype(joint_positions),2}}(undef, n + 1)
    Vi = zeros(eltype(joint_positions), 6, n + 1)
    Vdi = zeros(eltype(joint_positions), 6, n + 1)
    Vdi[:, 1] = vcat(zeros(eltype(joint_positions), 3), -gravity)
    AdTi[n+1] = adjoint_representation(transform_inv(link_frames[n+1]))
    Fi = copy(tip_wrench)
    joint_torques = zeros(eltype(joint_positions), n)

    for i = 1:n
        Mi *= link_frames[i]
        Ai[:, i] = adjoint_representation(transform_inv(Mi)) * screw_axes[:, i]
        AdTi[i] = adjoint_representation(
            matrix_exp6(vec_to_se3(Ai[:, i] * -joint_positions[i])) *
            transform_inv(link_frames[i]),
        )
        Vi[:, i+1] = AdTi[i] * Vi[:, i] + Ai[:, i] * joint_velocities[i]
        Vdi[:, i+1] =
            AdTi[i] * Vdi[:, i] +
            Ai[:, i] * joint_accelerations[i] +
            ad(Vi[:, i+1]) * Ai[:, i] * joint_velocities[i]
    end

    for i = n:-1:1
        Fi =
            AdTi[i+1]' * Fi + spatial_inertias[i] * Vdi[:, i+1] -
            ad(Vi[:, i+1])' * spatial_inertias[i] * Vi[:, i+1]
        joint_torques[i] = Fi' * Ai[:, i]
    end

    return joint_torques
end

"""
    mass_matrix(joint_positions, link_frames, spatial_inertias, screw_axes)

Computes the mass matrix of an open chain robot based on the given configuration.
Calls [`inverse_dynamics`](@ref) ``n`` times, each time passing a ``\\ddot{\\theta}``
vector with a single element equal to one and all other inputs set to zero. Each call
generates a single column of ``M(\\theta)``.

!!! info "What is the mass matrix?"
    The mass matrix ``M(\\theta)`` relates joint accelerations to joint torques in the
    absence of velocity-dependent and gravitational forces:
    ``\\tau = M(\\theta) \\ddot{\\theta}``. It captures how the robot's inertia is
    distributed across its joints at a given configuration. This function computes it by
    calling [`inverse_dynamics`](@ref) ``n`` times, once per joint.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}``.
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.

# Returns
The ``n×n`` mass matrix ``M(\\theta)``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> mass_matrix([0.1, 0.1, 0.1], link_frames, spatial_inertias, screw_axes)
3×3 Matrix{Float64}:
 22.5433      -0.307147  -0.00718426
 -0.307147     1.96851    0.432157
 -0.00718426   0.432157   0.191631
```
"""
function mass_matrix(
    joint_positions::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    M = zeros(n, n)

    for i = 1:n
        joint_accelerations = zeros(n)
        joint_accelerations[i] = 1
        M[:, i] = inverse_dynamics(
            joint_positions,
            zeros(n),
            joint_accelerations,
            zeros(3),
            zeros(6),
            link_frames,
            spatial_inertias,
            screw_axes,
        )
    end

    return M
end

"""
    velocity_quadratic_forces(joint_positions, joint_velocities, link_frames, spatial_inertias, screw_axes)

Computes the Coriolis and centripetal terms ``c(\\theta, \\dot{\\theta})`` in the inverse
dynamics of an open chain robot. Calls [`inverse_dynamics`](@ref) with `gravity = 0`,
`tip_wrench = 0`, and `joint_accelerations = 0`.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `joint_velocities`: the ``n``-vector of joint rates.
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}``.
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.

# Returns
The ``n``-vector of Coriolis and centripetal terms.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> velocity_quadratic_forces([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], link_frames, spatial_inertias, screw_axes)
3-element Vector{Float64}:
  0.26453118054501235
 -0.0550515682891655
 -0.006891320068248912
```
"""
function velocity_quadratic_forces(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    inverse_dynamics(
        joint_positions,
        joint_velocities,
        zeros(length(joint_positions)),
        zeros(3),
        zeros(6),
        link_frames,
        spatial_inertias,
        screw_axes,
    )
end

"""
    gravity_forces(joint_positions, gravity, link_frames, spatial_inertias, screw_axes)

Computes the joint forces/torques an open chain robot requires to overcome gravity at
its configuration. Calls [`inverse_dynamics`](@ref) with `joint_velocities = 0`,
`joint_accelerations = 0`, and `tip_wrench = 0`.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `gravity`: the gravity vector ``g`` (e.g., `[0, 0, -9.8]`).
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}``.
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.

# Returns
The ``n``-vector of joint gravity torques ``g(\\theta)``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> gravity_forces([0.1, 0.1, 0.1], [0, 0, -9.8], link_frames, spatial_inertias, screw_axes)
3-element Vector{Float64}:
  28.40331261821983
 -37.64094817177068
  -5.4415891999683605
```
"""
function gravity_forces(
    joint_positions::AbstractVector,
    gravity::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    inverse_dynamics(
        joint_positions,
        zeros(n),
        zeros(n),
        gravity,
        zeros(6),
        link_frames,
        spatial_inertias,
        screw_axes,
    )
end

"""
    end_effector_forces(joint_positions, tip_wrench, link_frames, spatial_inertias, screw_axes)

Computes the joint forces/torques an open chain robot requires only to create the end-effector force `tip_wrench`.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `tip_wrench`: the spatial force applied by the end-effector expressed in frame `{n+1}`.
- `link_frames`: the list of link frames `i` relative to `i-1` at the home position.
- `spatial_inertias`: the spatial inertia matrices `Gi` of the links.
- `screw_axes`: the screw axes `Si` of the joints in a space frame, in the format of a matrix with axes as the columns.

Returns the joint forces and torques required only to create the end-effector force `tip_wrench`.
This function calls inverse_dynamics with `gravity = 0`, `joint_velocities = 0`, and `joint_accelerations = 0`.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> import LinearAlgebra as LA

julia> joint_positions = [0.1, 0.1, 0.1]
3-element Vector{Float64}:
 0.1
 0.1
 0.1

julia> tip_wrench = [1, 1, 1, 1, 1, 1]
6-element Vector{Int64}:
 1
 1
 1
 1
 1
 1

julia> M01 = [1 0 0        0;
              0 1 0        0;
              0 0 1 0.089159;
              0 0 0        1]
4×4 Matrix{Float64}:
 1.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.089159
 0.0  0.0  0.0  1.0

julia> M12 = [ 0 0 1    0.28;
               0 1 0 0.13585;
              -1 0 0       0;
               0 0 0       1]
4×4 Matrix{Float64}:
  0.0  0.0  1.0  0.28
  0.0  1.0  0.0  0.13585
 -1.0  0.0  0.0  0.0
  0.0  0.0  0.0  1.0

julia> M23 = [1 0 0       0;
              0 1 0 -0.1197;
              0 0 1   0.395;
              0 0 0       1]
4×4 Matrix{Float64}:
 1.0  0.0  0.0   0.0
 0.0  1.0  0.0  -0.1197
 0.0  0.0  1.0   0.395
 0.0  0.0  0.0   1.0

julia> M34 = [1 0 0       0;
              0 1 0       0;
              0 0 1 0.14225;
              0 0 0       1]
4×4 Matrix{Float64}:
 1.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.14225
 0.0  0.0  0.0  1.0

julia> link_frames = [M01, M12, M23, M34]
4-element Vector{Matrix{Float64}}:
 [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.089159; 0.0 0.0 0.0 1.0]
 [0.0 0.0 1.0 0.28; 0.0 1.0 0.0 0.13585; -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0]
 [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 -0.1197; 0.0 0.0 1.0 0.395; 0.0 0.0 0.0 1.0]
 [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.14225; 0.0 0.0 0.0 1.0]

julia> G1 = LA.Diagonal([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
6×6 LinearAlgebra.Diagonal{Float64, Vector{Float64}}:
 0.010267   ⋅         ⋅        ⋅    ⋅    ⋅ 
  ⋅        0.010267   ⋅        ⋅    ⋅    ⋅ 
  ⋅         ⋅        0.00666   ⋅    ⋅    ⋅ 
  ⋅         ⋅         ⋅       3.7   ⋅    ⋅ 
  ⋅         ⋅         ⋅        ⋅   3.7   ⋅ 
  ⋅         ⋅         ⋅        ⋅    ⋅   3.7

julia> G2 = LA.Diagonal([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
6×6 LinearAlgebra.Diagonal{Float64, Vector{Float64}}:
 0.22689   ⋅        ⋅          ⋅      ⋅      ⋅ 
  ⋅       0.22689   ⋅          ⋅      ⋅      ⋅ 
  ⋅        ⋅       0.0151074   ⋅      ⋅      ⋅ 
  ⋅        ⋅        ⋅         8.393   ⋅      ⋅ 
  ⋅        ⋅        ⋅          ⋅     8.393   ⋅ 
  ⋅        ⋅        ⋅          ⋅      ⋅     8.393

julia> G3 = LA.Diagonal([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])
6×6 LinearAlgebra.Diagonal{Float64, Vector{Float64}}:
 0.0494433   ⋅          ⋅         ⋅      ⋅      ⋅ 
  ⋅         0.0494433   ⋅         ⋅      ⋅      ⋅ 
  ⋅          ⋅         0.004095   ⋅      ⋅      ⋅ 
  ⋅          ⋅          ⋅        2.275   ⋅      ⋅ 
  ⋅          ⋅          ⋅         ⋅     2.275   ⋅ 
  ⋅          ⋅          ⋅         ⋅      ⋅     2.275

julia> spatial_inertias = [G1, G2, G3]
3-element Vector{LinearAlgebra.Diagonal{Float64, Vector{Float64}}}:
 Diagonal([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
 Diagonal([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
 Diagonal([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])

julia> screw_axes = [ 1  0  1      0  1      0;
                 0  1  0 -0.089  0      0;
                 0  1  0 -0.089  0  0.425]'
6×3 adjoint(::Matrix{Float64}) with eltype Float64:
 1.0   0.0     0.0
 0.0   1.0     1.0
 1.0   0.0     0.0
 0.0  -0.089  -0.089
 1.0   0.0     0.0
 0.0   0.0     0.425

julia> end_effector_forces(joint_positions, tip_wrench, link_frames, spatial_inertias, screw_axes)
3-element Vector{Float64}:
 1.4095460782639782
 1.8577149723180628
 1.392409
```
"""
function end_effector_forces(
    joint_positions::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    inverse_dynamics(
        joint_positions,
        zeros(n),
        zeros(n),
        zeros(3),
        tip_wrench,
        link_frames,
        spatial_inertias,
        screw_axes,
    )
end

"""
    forward_dynamics(joint_positions, joint_velocities, joint_torques, gravity, tip_wrench, link_frames, spatial_inertias, screw_axes)

Computes forward dynamics in the space frame for an open chain robot.
Computes ``\\ddot{\\theta}`` by solving
``M(\\theta)\\ddot{\\theta} = \\tau - c(\\theta,\\dot{\\theta}) - g(\\theta) - J^T(\\theta) \\mathcal{F}_{\\text{tip}}``.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `joint_velocities`: the ``n``-vector of joint rates.
- `joint_torques`: the ``n``-vector of joint forces/torques.
- `gravity`: the gravity vector ``g`` (e.g., `[0, 0, -9.8]`).
- `tip_wrench`: the wrench ``\\mathcal{F}_{\\text{tip}}`` applied by the end-effector expressed in frame ``\\{n+1\\}``.
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}``.
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.

# Returns
The ``n``-vector of joint accelerations ``\\ddot{\\theta}``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> forward_dynamics([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], [0.5, 0.6, 0.7], [0, 0, -9.8], [1, 1, 1, 1, 1, 1], link_frames, spatial_inertias, screw_axes)
3-element Vector{Float64}:
  -0.9739290670855625
  25.584667840340547
 -32.91499212478147
```
"""
function forward_dynamics(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_torques::AbstractVector,
    gravity::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    LA.inv(mass_matrix(joint_positions, link_frames, spatial_inertias, screw_axes)) * (
        joint_torques - velocity_quadratic_forces(
            joint_positions,
            joint_velocities,
            link_frames,
            spatial_inertias,
            screw_axes,
        ) - gravity_forces(
            joint_positions,
            gravity,
            link_frames,
            spatial_inertias,
            screw_axes,
        ) - end_effector_forces(
            joint_positions,
            tip_wrench,
            link_frames,
            spatial_inertias,
            screw_axes,
        )
    )
end

"""
    euler_step(joint_positions, joint_velocities, joint_accelerations, timestep)

Compute the joint angles and velocities at the next timestep using first order Euler integration.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `joint_velocities`: the ``n``-vector of joint rates.
- `joint_accelerations`: the ``n``-vector of joint accelerations.
- `timestep`: the timestep delta t.

# Return
- `joint_positionsNext`: the vector of joint variables after `timestep` from first order Euler integration.
- `joint_velocitiesNext`: the vector of joint rates after `timestep` from first order Euler integration.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> euler_step([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], [2, 1.5, 1], 0.1)
([0.11000000000000001, 0.12000000000000001, 0.13], [0.30000000000000004, 0.35000000000000003, 0.4])
```
"""
function euler_step(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_accelerations::AbstractVector,
    timestep::Number,
)
    joint_positions + timestep * joint_velocities,
    joint_velocities + timestep * joint_accelerations
end

"""
    inverse_dynamics_trajectory(joint_position_traj, joint_velocity_traj, joint_acceleration_traj, gravity, tip_wrench_traj, link_frames, spatial_inertias, screw_axes)

Calculates the joint forces/torques required to move the serial chain along the given
trajectory using [`inverse_dynamics`](@ref) at each timestep.

# Arguments
- `joint_position_traj`: an ``N×n`` matrix of joint variables, where row ``i`` is the joint vector at timestep ``i``.
- `joint_velocity_traj`: an ``N×n`` matrix of joint velocities.
- `joint_acceleration_traj`: an ``N×n`` matrix of joint accelerations.
- `gravity`: the gravity vector ``g`` (e.g., `[0, 0, -9.8]`).
- `tip_wrench_traj`: an ``N×6`` matrix where row ``i`` is the spatial wrench applied by the end-effector at timestep ``i``.
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}``.
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.

# Returns
An ``N×n`` matrix of joint forces/torques, where row ``i`` is the joint force-torque at timestep ``i``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> traj_pos = [0.0 0.0 0.0; 0.5 0.5 0.5; 1.0 1.0 1.0];

julia> traj_vel = [0.0 0.0 0.0; 0.5 0.5 0.5; 1.0 1.0 1.0];

julia> traj_acc = [0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5];

julia> inverse_dynamics_trajectory(traj_pos, traj_vel, traj_acc, [0, 0, -9.8], ones(3, 6), link_frames, spatial_inertias, screw_axes)
3×3 adjoint(::Matrix{Float64}) with eltype Float64:
  23.7113  -35.2321  -3.87345
 105.336   -27.6441  -1.41169
 129.143   -17.3084   2.30991
```
"""
function inverse_dynamics_trajectory(
    joint_position_traj::AbstractMatrix,
    joint_velocity_traj::AbstractMatrix,
    joint_acceleration_traj::AbstractMatrix,
    gravity::AbstractVector,
    tip_wrench_traj::AbstractMatrix,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    joint_position_traj = joint_position_traj'
    joint_velocity_traj = joint_velocity_traj'
    joint_acceleration_traj = joint_acceleration_traj'
    tip_wrench_traj = tip_wrench_traj'
    joint_torque_traj = copy(joint_position_traj)

    for i in axes(joint_position_traj, 2)
        joint_torque_traj[:, i] = inverse_dynamics(
            joint_position_traj[:, i],
            joint_velocity_traj[:, i],
            joint_acceleration_traj[:, i],
            gravity,
            tip_wrench_traj[:, i],
            link_frames,
            spatial_inertias,
            screw_axes,
        )
    end

    joint_torque_traj'
end

"""
    forward_dynamics_trajectory(joint_positions, joint_velocities, joint_torque_traj, gravity, tip_wrench_traj, link_frames, spatial_inertias, screw_axes, timestep, integration_resolution)

Simulates the motion of a serial chain given an open-loop history of joint
forces/torques. Uses [`forward_dynamics`](@ref) with [`euler_step`](@ref) integration at
each timestep.

# Arguments
- `joint_positions`: the ``n``-vector of initial joint variables.
- `joint_velocities`: the ``n``-vector of initial joint velocities.
- `joint_torque_traj`: an ``N×n`` matrix where row ``i`` is the joint force/torque vector at timestep ``i``.
- `gravity`: the gravity vector ``g`` (e.g., `[0, 0, -9.8]`).
- `tip_wrench_traj`: an ``N×6`` matrix where row ``i`` is the spatial wrench applied by the end-effector at timestep ``i``.
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}``.
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.
- `timestep`: the timestep ``\\Delta t`` between consecutive trajectory points.
- `integration_resolution`: the number of Euler integration steps during each timestep ``\\Delta t``. Must be a positive integer ``\\geq 1``. Larger values result in slower simulations but less accumulation of integration error.

# Returns
- `joint_position_traj`: the ``N×n`` matrix of joint variables resulting from the specified joint forces/torques.
- `joint_velocity_traj`: the ``N×n`` matrix of joint velocities resulting from the specified joint forces/torques.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> joint_torque_traj = [3.63 -6.58 -5.57; 3.74 -5.55 -5.5; 4.31 -0.68 -5.19];

julia> thetamat, dthetamat = forward_dynamics_trajectory([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], joint_torque_traj, [0, 0, -9.8], ones(3, 6), link_frames, spatial_inertias, screw_axes, 0.1, 8);

julia> thetamat
3×3 adjoint(::Matrix{Float64}) with eltype Float64:
 0.1       0.1        0.1
 0.106431  0.2626    -0.226649
 0.10198   0.715813  -1.22522
```
"""
function forward_dynamics_trajectory(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_torque_traj::AbstractMatrix,
    gravity::AbstractVector,
    tip_wrench_traj::AbstractMatrix,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
    timestep::Number,
    integration_resolution::Number,
)
    joint_torque_traj = joint_torque_traj'
    tip_wrench_traj = tip_wrench_traj'
    joint_position_traj = copy(joint_torque_traj)
    joint_position_traj[:, 1] = joint_positions
    joint_velocity_traj = copy(joint_torque_traj)
    joint_velocity_traj[:, 1] = joint_velocities

    for i = first(axes(joint_torque_traj, 2)):(last(axes(joint_torque_traj, 2))-1)
        for j = 1:integration_resolution
            joint_accelerations = forward_dynamics(
                joint_positions,
                joint_velocities,
                joint_torque_traj[:, i],
                gravity,
                tip_wrench_traj[:, i],
                link_frames,
                spatial_inertias,
                screw_axes,
            )
            joint_positions, joint_velocities = euler_step(
                joint_positions,
                joint_velocities,
                joint_accelerations,
                1.0 * timestep / integration_resolution,
            )
        end

        joint_position_traj[:, i+1] = joint_positions
        joint_velocity_traj[:, i+1] = joint_velocities
    end

    joint_position_traj = joint_position_traj'
    joint_velocity_traj = joint_velocity_traj'
    return joint_position_traj, joint_velocity_traj
end

# """
# *** CHAPTER 9: TRAJECTORY GENERATION ***
# """

"""
    cubic_time_scaling(total_time, t)

Computes ``s(t) = 3(t/T_f)^2 - 2(t/T_f)^3`` for a cubic time scaling, satisfying ``s(0)=0``, ``s(T_f)=1``, ``\\dot{s}(0)=0``, ``\\dot{s}(T_f)=0``.

# Arguments
- `total_time`: the total time ``T_f`` of the motion in seconds.
- `t`: the current time ``t`` satisfying ``0 \\le t \\le T_f``.

# Returns
A scalar ``s \\in [0, 1]`` representing the fraction of motion completed.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> cubic_time_scaling(2, 0.6)
0.21600000000000003
```
"""
cubic_time_scaling(total_time::Number, t::Number) =
    3(t / total_time)^2 - 2(t / total_time)^3

"""
    quintic_time_scaling(total_time, t)

Computes ``s(t) = 10(t/T_f)^3 - 15(t/T_f)^4 + 6(t/T_f)^5`` for a quintic time scaling, satisfying ``s(0)=0``, ``s(T_f)=1``, ``\\dot{s}(0)=0``, ``\\dot{s}(T_f)=0``, ``\\ddot{s}(0)=0``, ``\\ddot{s}(T_f)=0``.

# Arguments
- `total_time`: the total time ``T_f`` of the motion in seconds.
- `t`: the current time ``t`` satisfying ``0 \\le t \\le T_f``.

# Returns
A scalar ``s \\in [0, 1]`` representing the fraction of motion completed.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> quintic_time_scaling(2, 0.6)
0.16308
```
"""
quintic_time_scaling(total_time::Number, t::Number) =
    10(t / total_time)^3 - 15(t / total_time)^4 + 6(t / total_time)^5

"""
    joint_trajectory(joint_position_start, joint_position_end, total_time, N, method)

Computes a straight-line trajectory in joint space with the specified time scaling.

# Arguments
- `joint_position_start`: the ``n``-vector of initial joint variables.
- `joint_position_end`: the ``n``-vector of final joint variables.
- `total_time`: total time of the motion in seconds from rest to rest.
- `N`: the number of points ``N \\geq 2`` in the discrete representation of the trajectory.
- `method`: the time-scaling method — use `3` for cubic ([`cubic_time_scaling`](@ref)) or `5` for quintic ([`quintic_time_scaling`](@ref)).

# Returns
An ``N×n`` matrix where each row is an ``n``-vector of joint variables. The first row is `joint_position_start` and the ``N``th row is `joint_position_end`. The elapsed time between each row is ``T_f / (N - 1)``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> joint_trajectory([0, 0], [1, 2], 2, 4, 3)
4×2 adjoint(::Matrix{Float64}) with eltype Float64:
 0.0       0.0
 0.259259  0.518519
 0.740741  1.48148
 1.0       2.0
```
"""
function joint_trajectory(
    joint_position_start::AbstractVector,
    joint_position_end::AbstractVector,
    total_time::Number,
    N::Integer,
    method::Integer,
)
    timegap = total_time / (N - 1)
    traj = zeros(length(joint_position_start), N)

    for i = 1:N
        if method == 3
            s = cubic_time_scaling(total_time, timegap * (i - 1))
        else
            s = quintic_time_scaling(total_time, timegap * (i - 1))
        end

        traj[:, i] = s * joint_position_end + (1 - s) * joint_position_start
    end

    traj'
end

"""
    screw_trajectory(transform_start, transform_end, total_time, N, method)

Computes a trajectory as a list of ``N`` SE(3) matrices corresponding to the screw
motion about a space screw axis. Unlike [`cartesian_trajectory`](@ref), the origin
follows a screw path rather than a straight line.

!!! info "What is a screw trajectory?"
    A screw trajectory moves the end-effector along a constant screw axis — like turning
    a corkscrew. The rotation and translation are coupled: the end-effector follows a
    helical path. Compare with [`cartesian_trajectory`](@ref), which decouples them.

# Arguments
- `transform_start`: the initial end-effector configuration (4×4 SE(3) matrix).
- `transform_end`: the final end-effector configuration (4×4 SE(3) matrix).
- `total_time`: total time of the motion in seconds from rest to rest.
- `N`: the number of points ``N \\geq 2`` in the discrete representation of the trajectory.
- `method`: the time-scaling method — use `3` for cubic ([`cubic_time_scaling`](@ref)) or `5` for quintic ([`quintic_time_scaling`](@ref)).

# Returns
A list of ``N`` SE(3) matrices separated in time by ``T_f / (N - 1)``. The first in the list is `transform_start` and the ``N``th is `transform_end`.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> Xstart = [1 0 0 1; 0 1 0 0; 0 0 1 1; 0 0 0 1];

julia> Xend = [0 0 1 0.1; 1 0 0 0; 0 1 0 4.1; 0 0 0 1];

julia> traj = screw_trajectory(Xstart, Xend, 5, 4, 3);

julia> length(traj)
4

julia> traj[1]
4×4 Matrix{Float64}:
 1.0  0.0  0.0  1.0
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  1.0
 0.0  0.0  0.0  1.0
```
"""
function screw_trajectory(
    transform_start::AbstractMatrix,
    transform_end::AbstractMatrix,
    total_time::Number,
    N::Integer,
    method::Integer,
)
    timegap = total_time / (N - 1)
    traj = Array{Array{Float64,2}}(undef, N)

    for i = 1:N
        if method == 3
            s = cubic_time_scaling(total_time, timegap * (i - 1))
        else
            s = quintic_time_scaling(total_time, timegap * (i - 1))
        end

        traj[i] =
            transform_start *
            matrix_exp6(matrix_log6(transform_inv(transform_start) * transform_end) * s)
    end

    return traj
end

"""
    cartesian_trajectory(transform_start, transform_end, total_time, N, method)

Computes a trajectory as a list of ``N`` SE(3) matrices where the origin of the
end-effector frame follows a straight line and the rotation follows a geodesic on
SO(3). Unlike [`screw_trajectory`](@ref), the position is decoupled from the rotation.

!!! info "Screw vs Cartesian trajectory"
    Unlike [`screw_trajectory`](@ref), a Cartesian trajectory decouples rotation and
    translation: the origin follows a straight line in space while the orientation
    interpolates independently along the shortest path (geodesic) on SO(3). This is
    usually more intuitive for pick-and-place tasks.

# Arguments
- `transform_start`: the initial end-effector configuration (4×4 SE(3) matrix).
- `transform_end`: the final end-effector configuration (4×4 SE(3) matrix).
- `total_time`: total time of the motion in seconds from rest to rest.
- `N`: the number of points ``N \\geq 2`` in the discrete representation of the trajectory.
- `method`: the time-scaling method — use `3` for cubic ([`cubic_time_scaling`](@ref)) or `5` for quintic ([`quintic_time_scaling`](@ref)).

# Returns
A list of ``N`` SE(3) matrices separated in time by ``T_f / (N - 1)``. The first in the list is `transform_start` and the ``N``th is `transform_end`.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> Xstart = [1 0 0 1; 0 1 0 0; 0 0 1 1; 0 0 0 1];

julia> Xend = [0 0 1 0.1; 1 0 0 0; 0 1 0 4.1; 0 0 0 1];

julia> traj = cartesian_trajectory(Xstart, Xend, 5, 4, 5);

julia> length(traj)
4

julia> traj[1]
4×4 Matrix{Float64}:
 1.0  0.0  0.0  1.0
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  1.0
 0.0  0.0  0.0  1.0
```
"""
function cartesian_trajectory(
    transform_start::AbstractMatrix,
    transform_end::AbstractMatrix,
    total_time::Number,
    N::Integer,
    method::Integer,
)
    timegap = total_time / (N - 1)
    traj = Array{Array{Float64,2}}(undef, N)

    Rstart, pstart = transform_to_rotation_position(transform_start)
    Rend, pend = transform_to_rotation_position(transform_end)

    for i = 1:N
        if method == 3
            s = cubic_time_scaling(total_time, timegap * (i - 1))
        else
            s = quintic_time_scaling(total_time, timegap * (i - 1))
        end

        traj[i] = vcat(
            hcat(
                Rstart * matrix_exp3(matrix_log3(Rstart' * Rend) * s),
                s * pend + (1 - s) * pstart,
            ),
            [0 0 0 1],
        )
    end

    return traj
end

# """
# *** CHAPTER 11: ROBOT CONTROL ***
# """

"""
    computed_torque(joint_positions, joint_velocities, error_integral, gravity, link_frames, spatial_inertias, screw_axes, desired_joint_positions, desired_joint_velocities, desired_joint_accelerations, Kp, Ki, Kd)

Computes the joint control torques at a particular time instant using the computed
torque control law:

``\\tau = \\widehat{M}(\\theta)(\\ddot{\\theta}_d + K_p e + K_i \\int e + K_d \\dot{e}) + \\widehat{h}(\\theta, \\dot{\\theta})``

where ``e = \\theta_d - \\theta``, ``\\dot{e} = \\dot{\\theta}_d - \\dot{\\theta}``,
``\\widehat{M}`` is the model of the robot's mass matrix, and ``\\widehat{h}`` is the
model of centripetal, Coriolis, and gravitational forces.

!!! info "What is computed torque control?"
    Computed torque control (also called inverse dynamics control) uses a full dynamic
    model of the robot to cancel out nonlinear effects (Coriolis, centripetal, gravity),
    then applies a PID-like feedback law on the linearised system. This achieves precise
    trajectory tracking even at high speeds, unlike simpler PD controllers that ignore
    the robot's dynamics.

# Arguments
- `joint_positions`: the ``n``-vector of current joint variables ``\\theta``.
- `joint_velocities`: the ``n``-vector of current joint rates ``\\dot{\\theta}``.
- `error_integral`: the ``n``-vector of the time-integral of joint errors ``\\int e \\, dt``.
- `gravity`: the gravity vector ``g`` (e.g., `[0, 0, -9.8]`).
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}``.
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.
- `desired_joint_positions`: the ``n``-vector of desired joint variables ``\\theta_d``.
- `desired_joint_velocities`: the ``n``-vector of desired joint rates ``\\dot{\\theta}_d``.
- `desired_joint_accelerations`: the ``n``-vector of desired joint accelerations ``\\ddot{\\theta}_d``.
- `Kp`: the proportional gain (scalar, applied as ``K_p I``).
- `Ki`: the integral gain (scalar).
- `Kd`: the derivative gain (scalar).

# Returns
The ``n``-vector of joint control torques.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> computed_torque([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], [0.2, 0.2, 0.2], [0, 0, -9.8], link_frames, spatial_inertias, screw_axes, [1.0, 1.0, 1.0], [2.0, 1.2, 2.0], [0.1, 0.1, 0.1], 1.3, 1.2, 1.1)
3-element Vector{Float64}:
 133.0052524649953
 -29.942233243760633
  -3.03276856161724
```
"""
function computed_torque(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    error_integral::AbstractVector,
    gravity::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
    desired_joint_positions::AbstractVector,
    desired_joint_velocities::AbstractVector,
    desired_joint_accelerations::AbstractVector,
    Kp::Number,
    Ki::Number,
    Kd::Number,
)
    e = desired_joint_positions - joint_positions
    mass_matrix(joint_positions, link_frames, spatial_inertias, screw_axes) * (
        Kp * e +
        Ki * (error_integral + e) +
        Kd * (desired_joint_velocities - joint_velocities)
    ) + inverse_dynamics(
        joint_positions,
        joint_velocities,
        desired_joint_accelerations,
        gravity,
        zeros(6),
        link_frames,
        spatial_inertias,
        screw_axes,
    )
end

"""
    simulate_control(joint_positions, joint_velocities, gravity, tip_wrench_traj, link_frames, spatial_inertias, screw_axes, desired_joint_position_traj, desired_joint_velocity_traj, desired_joint_acceleration_traj, estimated_gravity, estimated_link_frames, estimated_spatial_inertias, Kp, Ki, Kd, timestep, integration_resolution)

Simulates the [`computed_torque`](@ref) controller over a given desired trajectory using
[`forward_dynamics`](@ref) and numerical integration. The controller may use different
(possibly inaccurate) dynamics parameters
(`estimated_gravity`, `estimated_link_frames`, `estimated_spatial_inertias`) while the
actual forward dynamics simulation uses the true parameters (`gravity`, `link_frames`,
`spatial_inertias`).

!!! info "Why two sets of dynamics parameters?"
    In real robotics, the controller's model of the robot is never perfectly accurate.
    This function simulates that reality: the "actual" dynamics use the true parameters
    (`gravity`, `link_frames`, `spatial_inertias`), while the controller uses its own
    estimates (`gravity_estimate`, `link_frames_estimate`,
    `spatial_inertias_estimate`). The difference shows how robust the controller is to
    modelling errors.

# Arguments
- `joint_positions`: the ``n``-vector of initial joint variables.
- `joint_velocities`: the ``n``-vector of initial joint velocities.
- `gravity`: the actual gravity vector used in the forward dynamics simulation.
- `tip_wrench_traj`: an ``N×6`` matrix where row ``i`` is the spatial wrench applied by the end-effector at timestep ``i``.
- `link_frames`: the actual link frames (vector of ``n+1`` SE(3) matrices) used in the forward dynamics simulation.
- `spatial_inertias`: the actual spatial inertia matrices (vector of ``n`` 6×6 matrices) used in the forward dynamics simulation.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.
- `desired_joint_position_traj`: an ``N×n`` matrix of desired joint positions at each timestep.
- `desired_joint_velocity_traj`: an ``N×n`` matrix of desired joint velocities at each timestep.
- `desired_joint_acceleration_traj`: an ``N×n`` matrix of desired joint accelerations at each timestep.
- `estimated_gravity`: the gravity vector used by the controller (may differ from `gravity`).
- `estimated_link_frames`: the link frames used by the controller (may differ from `link_frames`).
- `estimated_spatial_inertias`: the spatial inertia matrices used by the controller (may differ from `spatial_inertias`).
- `Kp`: the proportional gain (scalar).
- `Ki`: the integral gain (scalar).
- `Kd`: the derivative gain (scalar).
- `timestep`: the timestep ``\\Delta t`` between trajectory reference points.
- `integration_resolution`: the number of Euler integration steps per timestep.

# Returns
- `joint_torque_traj`: an ``N×n`` matrix of applied joint torques at each timestep.
- `joint_position_traj`: an ``N×n`` matrix of actual joint positions at each timestep.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> desired_pos = [0.1 0.1 0.1; 0.2 0.2 0.2; 0.3 0.3 0.3];

julia> desired_vel = [0.1 0.1 0.1; 0.1 0.1 0.1; 0.1 0.1 0.1];

julia> taumat, thetamat = simulate_control([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], [0, 0, -9.8], ones(3, 6), link_frames, spatial_inertias, screw_axes, desired_pos, desired_vel, zeros(3, 3), [0, 0, -9.8], link_frames, spatial_inertias, 20, 10, 18, 0.1, 4);

julia> taumat
3×3 adjoint(::Matrix{Float64}) with eltype Float64:
 29.2466  -42.7951  -6.91623
 93.5113  -24.4938  -0.00376585
 45.3612  -39.6324  -5.62033
```
"""
function simulate_control(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    gravity::AbstractVector,
    tip_wrench_traj::AbstractMatrix,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
    desired_joint_position_traj::AbstractMatrix,
    desired_joint_velocity_traj::AbstractMatrix,
    desired_joint_acceleration_traj::AbstractMatrix,
    estimated_gravity::AbstractVector,
    estimated_link_frames::AbstractVector,
    estimated_spatial_inertias::AbstractVector,
    Kp::Number,
    Ki::Number,
    Kd::Number,
    timestep::Number,
    integration_resolution::Number,
)
    tip_wrench_traj = tip_wrench_traj'
    desired_joint_position_traj = desired_joint_position_traj'
    desired_joint_velocity_traj = desired_joint_velocity_traj'
    desired_joint_acceleration_traj = desired_joint_acceleration_traj'

    m, n = size(desired_joint_position_traj)

    current_positions = copy(joint_positions)
    current_velocities = copy(joint_velocities)

    error_integral = reshape(zeros(m, 1), (m,))
    joint_torque_traj = zeros(size(desired_joint_position_traj))
    joint_position_traj = zeros(size(desired_joint_position_traj))

    for i = 1:n
        joint_torques = computed_torque(
            current_positions,
            current_velocities,
            error_integral,
            estimated_gravity,
            estimated_link_frames,
            estimated_spatial_inertias,
            screw_axes,
            desired_joint_position_traj[:, i],
            desired_joint_velocity_traj[:, i],
            desired_joint_acceleration_traj[:, i],
            Kp,
            Ki,
            Kd,
        )

        for j = 1:integration_resolution
            joint_accelerations = forward_dynamics(
                current_positions,
                current_velocities,
                joint_torques,
                gravity,
                tip_wrench_traj[:, i],
                link_frames,
                spatial_inertias,
                screw_axes,
            )
            current_positions, current_velocities = euler_step(
                current_positions,
                current_velocities,
                joint_accelerations,
                timestep / integration_resolution,
            )
        end

        joint_torque_traj[:, i] = joint_torques
        joint_position_traj[:, i] = current_positions

        error_integral += timestep * (desired_joint_position_traj[:, i] - current_positions)
    end

    joint_torque_traj', joint_position_traj'
end

include("robot.jl")

end # module

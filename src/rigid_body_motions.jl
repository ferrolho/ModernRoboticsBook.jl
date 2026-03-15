export vec_to_so3,
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
    is_se3

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
3×3 StaticArraysCore.SMatrix{3, 3, Int64, 9} with indices SOneTo(3)×SOneTo(3):
  0  -3   2
  3   0  -1
 -2   1   0
```
"""
function vec_to_so3(ω::AbstractVector)
    @SMatrix [
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
3-element StaticArraysCore.SVector{3, Int64} with indices SOneTo(3):
 1
 2
 3
```
"""
function so3_to_vec(so3mat::AbstractMatrix)
    SA[so3mat[3, 2], so3mat[1, 3], so3mat[2, 1]]
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
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -0.694921   0.713521  0.0892929
 -0.192007  -0.303785  0.933192
  0.692978   0.63135   0.348107
```
"""
function matrix_exp3(so3mat::AbstractMatrix)
    ωθ = so3_to_vec(so3mat)
    if near_zero(LA.norm(ωθ))
        return SMatrix{3,3,Float64}(LA.I)
    else
        θ = axis_angle3(ωθ)[2]
        ωmat = so3mat / θ
        return SMatrix{3,3,Float64}(LA.I) + sin(θ) * ωmat + (1 - cos(θ)) * ωmat * ωmat
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
        return @SMatrix zeros(3, 3)
    elseif acosinput <= -1
        if !near_zero(1 + R[3, 3])
            ω = (1 / √(2 * (1 + R[3, 3]))) * SA[R[1, 3], R[2, 3], 1+R[3, 3]]
        elseif !near_zero(1 + R[2, 2])
            ω = (1 / √(2 * (1 + R[2, 2]))) * SA[R[1, 2], 1+R[2, 2], R[3, 2]]
        else
            ω = (1 / √(2 * (1 + R[1, 1]))) * SA[1+R[1, 1], R[2, 1], R[3, 1]]
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
4×4 StaticArraysCore.SMatrix{4, 4, Float64, 16} with indices SOneTo(4)×SOneTo(4):
 1.0  0.0   0.0  1.0
 0.0  0.0  -1.0  2.0
 0.0  1.0   0.0  5.0
 0.0  0.0   0.0  1.0
```
"""
function rotation_position_to_transform(R::AbstractMatrix, p::AbstractVector)
    @SMatrix [
        R[1, 1] R[1, 2] R[1, 3] p[1]
        R[2, 1] R[2, 2] R[2, 3] p[2]
        R[3, 1] R[3, 2] R[3, 3] p[3]
        0.0 0.0 0.0 1.0
    ]
end

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
4×4 StaticArraysCore.SMatrix{4, 4, Float64, 16} with indices SOneTo(4)×SOneTo(4):
 1.0   0.0  0.0   0.0
 0.0   0.0  1.0  -3.0
 0.0  -1.0  0.0   0.0
 0.0   0.0  0.0   1.0
```
"""
function transform_inv(T::AbstractMatrix)
    R = SMatrix{3,3}(@view T[1:3, 1:3])
    p = SA[T[1, 4], T[2, 4], T[3, 4]]
    Rt = R'
    Rtp = Rt * p
    @SMatrix [
        Rt[1, 1] Rt[1, 2] Rt[1, 3] -Rtp[1]
        Rt[2, 1] Rt[2, 2] Rt[2, 3] -Rtp[2]
        Rt[3, 1] Rt[3, 2] Rt[3, 3] -Rtp[3]
        0.0 0.0 0.0 1.0
    ]
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
4×4 StaticArraysCore.SMatrix{4, 4, Float64, 16} with indices SOneTo(4)×SOneTo(4):
  0.0  -3.0   2.0  4.0
  3.0   0.0  -1.0  5.0
 -2.0   1.0   0.0  6.0
  0.0   0.0   0.0  0.0
```
"""
function vec_to_se3(V::AbstractVector)
    ω1, ω2, ω3, v1, v2, v3 = V
    @SMatrix [
        0.0 -ω3 ω2 v1
        ω3 0.0 -ω1 v2
        -ω2 ω1 0.0 v3
        0.0 0.0 0.0 0.0
    ]
end

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
6-element StaticArraysCore.SVector{6, Int64} with indices SOneTo(6):
 1
 2
 3
 4
 5
 6
```
"""
se3_to_vec(se3mat::AbstractMatrix) =
    SA[se3mat[3, 2], se3mat[1, 3], se3mat[2, 1], se3mat[1, 4], se3mat[2, 4], se3mat[3, 4]]

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
6×6 StaticArraysCore.SMatrix{6, 6, Int64, 36} with indices SOneTo(6)×SOneTo(6):
 1  0   0  0  0   0
 0  0  -1  0  0   0
 0  1   0  0  0   0
 0  0   3  1  0   0
 3  0   0  0  0  -1
 0  0   0  0  1   0
```
"""
function adjoint_representation(T::AbstractMatrix)
    R = SMatrix{3,3}(@view T[1:3, 1:3])
    p = SA[T[1, 4], T[2, 4], T[3, 4]]
    pR = vec_to_so3(p) * R
    z = zero(eltype(T))
    @SMatrix [
        R[1, 1] R[1, 2] R[1, 3] z z z
        R[2, 1] R[2, 2] R[2, 3] z z z
        R[3, 1] R[3, 2] R[3, 3] z z z
        pR[1, 1] pR[1, 2] pR[1, 3] R[1, 1] R[1, 2] R[1, 3]
        pR[2, 1] pR[2, 2] pR[2, 3] R[2, 1] R[2, 2] R[2, 3]
        pR[3, 1] pR[3, 2] pR[3, 3] R[3, 1] R[3, 2] R[3, 3]
    ]
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
6-element StaticArraysCore.SVector{6, Int64} with indices SOneTo(6):
  0
  0
  1
  0
 -3
  2
```
"""
function screw_to_axis(q::AbstractVector, s::AbstractVector, h::Number)
    v = LA.cross(q, s) + h * s
    SA[s[1], s[2], s[3], v[1], v[2], v[3]]
end

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
4×4 StaticArraysCore.SMatrix{4, 4, Float64, 16} with indices SOneTo(4)×SOneTo(4):
 1.0  0.0         0.0        0.0
 0.0  6.7949e-9  -1.0        1.01923e-8
 0.0  1.0         6.7949e-9  3.0
 0.0  0.0         0.0        1.0
```
"""
function matrix_exp6(se3mat::AbstractMatrix)
    omgtheta = SMatrix{3,3}(
        se3mat[1, 1],
        se3mat[2, 1],
        se3mat[3, 1],
        se3mat[1, 2],
        se3mat[2, 2],
        se3mat[3, 2],
        se3mat[1, 3],
        se3mat[2, 3],
        se3mat[3, 3],
    )
    ωθ = so3_to_vec(omgtheta)
    v = SA[se3mat[1, 4], se3mat[2, 4], se3mat[3, 4]]
    I3 = SMatrix{3,3,Float64}(LA.I)
    if near_zero(LA.norm(ωθ))
        R = I3
        p = v
    else
        θ = axis_angle3(ωθ)[2]
        ωmat = omgtheta / θ
        R = I3 + sin(θ) * ωmat + (1 - cos(θ)) * ωmat * ωmat
        p = (I3 * θ + (1 - cos(θ)) * ωmat + (θ - sin(θ)) * ωmat * ωmat) * v / θ
    end
    @SMatrix [
        R[1, 1] R[1, 2] R[1, 3] p[1]
        R[2, 1] R[2, 2] R[2, 3] p[2]
        R[3, 1] R[3, 2] R[3, 3] p[3]
        0.0 0.0 0.0 1.0
    ]
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
4×4 StaticArraysCore.SMatrix{4, 4, Float64, 16} with indices SOneTo(4)×SOneTo(4):
 0.0  0.0      0.0     0.0
 0.0  0.0     -1.5708  2.35619
 0.0  1.5708   0.0     2.35619
 0.0  0.0      0.0     0.0
```
"""
function matrix_log6(T::AbstractMatrix)
    R = SMatrix{3,3}(@view T[1:3, 1:3])
    ωmat = matrix_log3(R)
    I3 = SMatrix{3,3,Float64}(LA.I)
    v = SA[T[1, 4], T[2, 4], T[3, 4]]
    if iszero(ωmat)
        @SMatrix [
            0.0 0.0 0.0 v[1]
            0.0 0.0 0.0 v[2]
            0.0 0.0 0.0 v[3]
            0.0 0.0 0.0 0.0
        ]
    else
        θ = acos((LA.tr(R) - 1) / 2)
        Ginv = I3 - ωmat / 2 + (1 / θ - 1 / tan(θ / 2) / 2) * ωmat * ωmat / θ
        p = Ginv * v
        @SMatrix [
            ωmat[1, 1] ωmat[1, 2] ωmat[1, 3] p[1]
            ωmat[2, 1] ωmat[2, 2] ωmat[2, 3] p[2]
            ωmat[3, 1] ωmat[3, 2] ωmat[3, 3] p[3]
            0.0 0.0 0.0 0.0
        ]
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
4×4 StaticArraysCore.SMatrix{4, 4, Float64, 16} with indices SOneTo(4)×SOneTo(4):
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

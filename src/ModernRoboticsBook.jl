module ModernRoboticsBook

import LinearAlgebra as LA

export NearZero,
    Normalize,
    RotInv,
    VecToso3,
    so3ToVec,
    AxisAng3,
    MatrixExp3,
    MatrixLog3,
    RpToTrans,
    TransToRp,
    TransInv,
    VecTose3,
    se3ToVec,
    Adjoint,
    ScrewToAxis,
    AxisAng6,
    MatrixExp6,
    MatrixLog6,
    ProjectToSO3,
    ProjectToSE3,
    DistanceToSO3,
    DistanceToSE3,
    TestIfSO3,
    TestIfSE3,
    FKinBody,
    FKinSpace,
    JacobianBody,
    JacobianSpace,
    IKinBody,
    IKinSpace,
    ad,
    InverseDynamics,
    MassMatrix,
    VelQuadraticForces,
    GravityForces,
    EndEffectorForces,
    ForwardDynamics,
    EulerStep,
    InverseDynamicsTrajectory,
    ForwardDynamicsTrajectory,
    CubicTimeScaling,
    QuinticTimeScaling,
    JointTrajectory,
    ScrewTrajectory,
    CartesianTrajectory,
    ComputedTorque,
    SimulateControl

# """
# *** BASIC HELPER FUNCTIONS ***
# """

"""
    NearZero(z)

Determines whether a scalar is small enough to be treated as zero.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> NearZero(-1e-7)
true
```
"""
NearZero(z::Number) = abs(z) < 1e-6

"""
    Normalize(V)

Normalizes a vector.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> Normalize([1, 2, 3])
3-element Vector{Float64}:
 0.2672612419124244
 0.5345224838248488
 0.8017837257372732
```
"""
Normalize(V::AbstractVector) = V / LA.norm(V)

# """
# *** CHAPTER 3: RIGID-BODY MOTIONS ***
# """

"""
    RotInv(R)

Inverts a rotation matrix.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> RotInv([0 0 1; 1 0 0; 0 1 0])
3×3 adjoint(::Matrix{Int64}) with eltype Int64:
 0  1  0
 0  0  1
 1  0  0
```
"""
RotInv(R::AbstractMatrix) = R'

"""
    VecToso3(ω)

Converts a 3-vector to an so(3) representation.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> VecToso3([1, 2, 3])
3×3 Matrix{Int64}:
  0  -3   2
  3   0  -1
 -2   1   0
```
"""
function VecToso3(ω::AbstractVector)
    [
        0 -ω[3] ω[2]
        ω[3] 0 -ω[1]
        -ω[2] ω[1] 0
    ]
end

"""
    so3ToVec(so3mat)

Converts an so(3) representation to a 3-vector.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> so3ToVec([0 -3 2; 3 0 -1; -2 1 0])
3-element Vector{Int64}:
 1
 2
 3
```
"""
function so3ToVec(so3mat::AbstractMatrix)
    [so3mat[3, 2], so3mat[1, 3], so3mat[2, 1]]
end

"""
    AxisAng3(expc3)

Converts a 3-vector of exponential coordinates for rotation into axis-angle form.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> AxisAng3([1, 2, 3])
([0.2672612419124244, 0.5345224838248488, 0.8017837257372732], 3.7416573867739413)
```
"""
AxisAng3(expc3::AbstractVector) = Normalize(expc3), LA.norm(expc3)

"""
    MatrixExp3(so3mat)

Computes the matrix exponential of a matrix in so(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> MatrixExp3([0 -3 2; 3 0 -1; -2 1 0])
3×3 Matrix{Float64}:
 -0.694921   0.713521  0.0892929
 -0.192007  -0.303785  0.933192 
  0.692978   0.63135   0.348107 
```
"""
function MatrixExp3(so3mat::AbstractMatrix)
    ωθ = so3ToVec(so3mat)
    if NearZero(LA.norm(ωθ))
        return Matrix{Float64}(LA.I, 3, 3)
    else
        θ = AxisAng3(ωθ)[2]
        ωmat = so3mat / θ
        return LA.I + sin(θ) * ωmat + (1 - cos(θ)) * ωmat * ωmat
    end
end

"""
    MatrixLog3(R)

Computes the matrix logarithm of a rotation matrix.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> MatrixLog3([0 0 1; 1 0 0; 0 1 0])
3×3 Matrix{Float64}:
  0.0     -1.2092   1.2092
  1.2092   0.0     -1.2092
 -1.2092   1.2092   0.0   
```
"""
function MatrixLog3(R::AbstractMatrix)
    acosinput = (LA.tr(R) - 1) / 2
    if acosinput >= 1
        return zeros(3, 3)
    elseif acosinput <= -1
        if !NearZero(1 + R[3, 3])
            ω = (1 / √(2 * (1 + R[3, 3]))) * [R[1, 3], R[2, 3], 1 + R[3, 3]]
        elseif !NearZero(1 + R[2, 2])
            ω = (1 / √(2 * (1 + R[2, 2]))) * [R[1, 2], 1 + R[2, 2], R[3, 2]]
        else
            ω = (1 / √(2 * (1 + R[1, 1]))) * [1 + R[1, 1], R[2, 1], R[3, 1]]
        end
        return VecToso3(π * ω)
    else
        θ = acos(acosinput)
        return θ / 2 / sin(θ) * (R - R')
    end
end

"""
    RpToTrans(R, p)

Converts a rotation matrix and a position vector into homogeneous transformation matrix.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> RpToTrans([1 0 0; 0 0 -1; 0 1 0], [1, 2, 5])
4×4 Matrix{Int64}:
 1  0   0  1
 0  0  -1  2
 0  1   0  5
 0  0   0  1
```
"""
RpToTrans(R::AbstractMatrix, p::AbstractVector) = vcat(hcat(R, p), [0 0 0 1])

"""
    TransToRp(T)

Converts a homogeneous transformation matrix into a rotation matrix and position vector.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> TransToRp([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
([1 0 0; 0 0 -1; 0 1 0], [0, 0, 3])
```
"""
TransToRp(T::AbstractMatrix) = T[1:3, 1:3], T[1:3, 4]

"""
    TransInv(T)

Inverts a homogeneous transformation matrix.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> TransInv([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
4×4 Matrix{Int64}:
 1   0  0   0
 0   0  1  -3
 0  -1  0   0
 0   0  0   1
```
"""
function TransInv(T::AbstractMatrix)
    R, p = TransToRp(T)
    vcat(hcat(R', -R' * p), [0 0 0 1])
end

"""
    VecTose3(V)

Converts a spatial velocity vector into a 4x4 matrix in se3.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> VecTose3([1, 2, 3, 4, 5, 6])
4×4 Matrix{Float64}:
  0.0  -3.0   2.0  4.0
  3.0   0.0  -1.0  5.0
 -2.0   1.0   0.0  6.0
  0.0   0.0   0.0  0.0
```
"""
VecTose3(V::AbstractVector) = vcat(hcat(VecToso3(V[1:3]), V[4:6]), zeros(1, 4))

"""
    se3ToVec(se3mat)

Converts an se3 matrix into a spatial velocity vector.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> se3ToVec([0 -3 2 4; 3 0 -1 5; -2 1 0 6; 0 0 0 0])
6-element Vector{Int64}:
 1
 2
 3
 4
 5
 6
```
"""
se3ToVec(se3mat::AbstractMatrix) =
    vcat([se3mat[3, 2], se3mat[1, 3], se3mat[2, 1]], se3mat[1:3, 4])

"""
    Adjoint(T)

Computes the adjoint representation of a homogeneous transformation matrix.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> Adjoint([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
6×6 Matrix{Float64}:
 1.0  0.0   0.0  0.0  0.0   0.0
 0.0  0.0  -1.0  0.0  0.0   0.0
 0.0  1.0   0.0  0.0  0.0   0.0
 0.0  0.0   3.0  1.0  0.0   0.0
 3.0  0.0   0.0  0.0  0.0  -1.0
 0.0  0.0   0.0  0.0  1.0   0.0
```
"""
function Adjoint(T::AbstractMatrix)
    R, p = TransToRp(T)
    vcat(hcat(R, zeros(3, 3)), hcat(VecToso3(p) * R, R))
end

"""
    ScrewToAxis(q, s, h)

Takes a parametric description of a screw axis and converts it to a normalized screw axis.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> ScrewToAxis([3; 0; 0], [0; 0; 1], 2)
6-element Vector{Int64}:
  0
  0
  1
  0
 -3
  2
```
"""
ScrewToAxis(q::AbstractVector, s::AbstractVector, h::Number) =
    vcat(s, LA.cross(q, s) + h * s)

"""
    AxisAng6(expc6)

Converts a 6-vector of exponential coordinates into screw axis-angle form.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> AxisAng6([1, 0, 0, 1, 2, 3])
([1.0, 0.0, 0.0, 1.0, 2.0, 3.0], 1.0)
```
"""
function AxisAng6(expc6::AbstractVector)
    θ = LA.norm(expc6[1:3])
    if NearZero(θ)
        θ = LA.norm(expc6[4:6])
    end
    expc6 / θ, θ
end

"""
    MatrixExp6(se3mat)

Computes the matrix exponential of an se3 representation of exponential coordinates.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> MatrixExp6([0 0 0 0; 0 0 -1.57079632 2.35619449; 0 1.57079632 0 2.35619449; 0 0 0 0])
4×4 Matrix{Float64}:
 1.0  0.0         0.0        0.0       
 0.0  6.7949e-9  -1.0        1.01923e-8
 0.0  1.0         6.7949e-9  3.0       
 0.0  0.0         0.0        1.0       
```
"""
function MatrixExp6(se3mat::AbstractMatrix)
    ωθ = so3ToVec(se3mat[1:3, 1:3])
    if NearZero(LA.norm(ωθ))
        return vcat(hcat(Matrix{Float64}(LA.I, 3, 3), se3mat[1:3, 4]), [0 0 0 1])
    else
        θ = AxisAng3(ωθ)[2]
        ωmat = se3mat[1:3, 1:3] / θ
        return vcat(
            hcat(
                MatrixExp3(se3mat[1:3, 1:3]),
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
    MatrixLog6(T)

Computes the matrix logarithm of a homogeneous transformation matrix.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> MatrixLog6([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
4×4 Matrix{Float64}:
 0.0  0.0      0.0     0.0    
 0.0  0.0     -1.5708  2.35619
 0.0  1.5708   0.0     2.35619
 0.0  0.0      0.0     0.0    
```
"""
function MatrixLog6(T::AbstractMatrix)
    R, p = TransToRp(T)
    ωmat = MatrixLog3(R)
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
    ProjectToSO3(mat)

Returns a projection of mat into SO(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> ProjectToSO3([0.675 0.150  0.720; 0.370 0.771 -0.511; -0.630 0.619  0.472])
3×3 Matrix{Float64}:
  0.679011  0.148945   0.718859
  0.373207  0.773196  -0.512723
 -0.632187  0.616428   0.469421
```
"""
function ProjectToSO3(mat::AbstractMatrix)
    F = LA.svd(mat)
    R = F.U * F.Vt
    if LA.det(R) < 0
        R = F.U * LA.Diagonal([1, 1, -1]) * F.Vt
    end
    return R
end

"""
    ProjectToSE3(mat)

Returns a projection of mat into SE(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> ProjectToSE3([0.675 0.150 0.720 1.2; 0.370 0.771 -0.511 5.4; -0.630 0.619 0.472 3.6; 0.003 0.002 0.010 0.9])
4×4 Matrix{Float64}:
  0.679011  0.148945   0.718859  1.2
  0.373207  0.773196  -0.512723  5.4
 -0.632187  0.616428   0.469421  3.6
  0.0       0.0        0.0       1.0
```
"""
ProjectToSE3(mat::AbstractMatrix) = RpToTrans(ProjectToSO3(mat[1:3, 1:3]), mat[1:3, 4])

"""
    DistanceToSO3(mat)

Returns the Frobenius norm to describe the distance of mat from the SO(3) manifold.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> DistanceToSO3([1.0 0.0 0.0; 0.0 0.1 -0.95; 0.0 1.0 0.1])
0.08835298523536149
```
"""
DistanceToSO3(mat::AbstractMatrix) = LA.det(mat) > 0 ? LA.norm(mat'mat - LA.I) : 1e+9

"""
    DistanceToSE3(mat)

Returns the Frobenius norm to describe the distance of mat from the SE(3) manifold.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> DistanceToSE3([1.0 0.0 0.0 1.2; 0.0 0.1 -0.95 1.5; 0.0 1.0 0.1 -0.9; 0.0 0.0 0.1 0.98])
0.13493053768513638
```
"""
function DistanceToSE3(mat::AbstractMatrix)
    matR = mat[1:3, 1:3]
    if LA.det(matR) > 0
        LA.norm(hcat(vcat(matR'matR, zeros(1, 3)), mat[4, :]) - LA.I)
    else
        1e+9
    end
end

"""
    TestIfSO3(mat)

Returns true if mat is close to or on the manifold SO(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> TestIfSO3([1.0 0.0 0.0; 0.0 0.1 -0.95; 0.0 1.0 0.1])
false
```
"""
TestIfSO3(mat::AbstractMatrix) = abs(DistanceToSO3(mat)) < 1e-3

"""
    TestIfSE3(mat)

Returns true if mat is close to or on the manifold SE(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> TestIfSE3([1.0 0.0 0.0 1.2; 0.0 0.1 -0.95 1.5; 0.0 1.0 0.1 -0.9; 0.0 0.0 0.1 0.98])
false
```
"""
TestIfSE3(mat::AbstractMatrix) = abs(DistanceToSE3(mat)) < 1e-3

# """
# *** CHAPTER 4: FORWARD KINEMATICS ***
# """

"""
    FKinBody(home_config, body_screw_axes, joint_positions)

Computes forward kinematics in the body frame for an open chain robot.

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

julia> FKinBody(home_config, body_screw_axes, joint_positions)
4×4 Matrix{Float64}:
 -1.14424e-17  1.0           0.0  -5.0
  1.0          1.14424e-17   0.0   4.0
  0.0          0.0          -1.0   1.68584
  0.0          0.0           0.0   1.0
```
"""
function FKinBody(
    home_config::AbstractMatrix,
    body_screw_axes::AbstractMatrix,
    joint_positions::AbstractVector,
)
    for i = 1:length(joint_positions)
        home_config *= MatrixExp6(VecTose3(body_screw_axes[:, i] * joint_positions[i]))
    end
    home_config
end

"""
    FKinSpace(home_config, screw_axes, joint_positions)

Computes forward kinematics in the space frame for an open chain robot.

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

julia> FKinSpace(home_config, screw_axes, joint_positions)
4×4 Matrix{Float64}:
 -1.14424e-17  1.0           0.0  -5.0
  1.0          1.14424e-17   0.0   4.0
  0.0          0.0          -1.0   1.68584
  0.0          0.0           0.0   1.0
```
"""
function FKinSpace(
    home_config::AbstractMatrix,
    screw_axes::AbstractMatrix,
    joint_positions::AbstractVector,
)
    for i = length(joint_positions):-1:1
        home_config =
            MatrixExp6(VecTose3(screw_axes[:, i] * joint_positions[i])) * home_config
    end
    home_config
end

# """
# *** CHAPTER 5: VELOCITY KINEMATICS AND STATICS ***
# """

"""
    JacobianBody(body_screw_axes, joint_positions)

Computes the body Jacobian for an open chain robot.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> body_screw_axes = [0 0 1   0 0.2 0.2;
                1 0 0   2   0   3;
                0 1 0   0   2   1;
                1 0 0 0.2 0.3 0.4]';

julia> joint_positions = [0.2, 1.1, 0.1, 1.2];

julia> JacobianBody(body_screw_axes, joint_positions)
6×4 Matrix{Float64}:
 -0.0452841  0.995004    0.0       1.0
  0.743593   0.0930486   0.362358  0.0
 -0.667097   0.0361754  -0.932039  0.0
  2.32586    1.66809     0.564108  0.2
 -1.44321    2.94561     1.43307   0.3
 -2.0664     1.82882    -1.58869   0.4
```
"""
function JacobianBody(body_screw_axes::AbstractMatrix, joint_positions::AbstractVector)
    T = LA.I
    Jb = copy(body_screw_axes)
    for i = length(joint_positions)-1:-1:1
        T *= MatrixExp6(VecTose3(body_screw_axes[:, i+1] * -joint_positions[i+1]))
        Jb[:, i] = Adjoint(T) * body_screw_axes[:, i]
    end
    Jb
end

"""
    JacobianSpace(screw_axes, joint_positions)

Computes the space Jacobian for an open chain robot.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> screw_axes = [0 0 1   0 0.2 0.2;
                1 0 0   2   0   3;
                0 1 0   0   2   1;
                1 0 0 0.2 0.3 0.4]';

julia> joint_positions = [0.2, 1.1, 0.1, 1.2];

julia> JacobianSpace(screw_axes, joint_positions)
6×4 Matrix{Float64}:
 0.0  0.980067  -0.0901156   0.957494 
 0.0  0.198669   0.444554    0.284876 
 1.0  0.0        0.891207   -0.0452841
 0.0  1.95219   -2.21635    -0.511615 
 0.2  0.436541  -2.43713     2.77536  
 0.2  2.96027    3.23573     2.22512  
```
"""
function JacobianSpace(screw_axes::AbstractMatrix, joint_positions::AbstractVector)
    T = LA.I
    Js = copy(screw_axes)
    for i = 2:length(joint_positions)
        T *= MatrixExp6(VecTose3(screw_axes[:, i-1] * joint_positions[i-1]))
        Js[:, i] = Adjoint(T) * screw_axes[:, i]
    end
    Js
end

# """
# *** CHAPTER 6: INVERSE KINEMATICS ***
# """

"""
    IKinBody(body_screw_axes, home_config, target_config, initial_guess, angular_tolerance, linear_tolerance)

Computes inverse kinematics in the body frame for an open chain robot.

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

julia> IKinBody(body_screw_axes, home_config, target_config, initial_guess, angular_tolerance, linear_tolerance)
([1.5707381937148923, 2.999666997382942, 3.141539129217613], true)
```
"""
function IKinBody(
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
    Vb = se3ToVec(
        MatrixLog6(
            TransInv(FKinBody(home_config, body_screw_axes, joint_positions)) *
            target_config,
        ),
    )
    err = LA.norm(Vb[1:3]) > angular_tolerance || LA.norm(Vb[4:6]) > linear_tolerance
    while err && i < maxiterations
        joint_positions += LA.pinv(JacobianBody(body_screw_axes, joint_positions)) * Vb
        i += 1
        Vb = se3ToVec(
            MatrixLog6(
                TransInv(FKinBody(home_config, body_screw_axes, joint_positions)) *
                target_config,
            ),
        )
        err = LA.norm(Vb[1:3]) > angular_tolerance || LA.norm(Vb[4:6]) > linear_tolerance
    end
    return joint_positions, !err
end

"""
    IKinSpace(screw_axes, home_config, target_config, initial_guess, angular_tolerance, linear_tolerance)

Computes inverse kinematics in the space frame for an open chain robot.

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

julia> IKinSpace(screw_axes, home_config, target_config, initial_guess, angular_tolerance, linear_tolerance)
([1.5707378296567203, 2.999663844672524, 3.141534199856583], true)
```
"""
function IKinSpace(
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
    Tsb = FKinSpace(home_config, screw_axes, joint_positions)
    Vs = Adjoint(Tsb) * se3ToVec(MatrixLog6(TransInv(Tsb) * target_config))
    err = LA.norm(Vs[1:3]) > angular_tolerance || LA.norm(Vs[4:6]) > linear_tolerance
    while err && i < maxiterations
        joint_positions += LA.pinv(JacobianSpace(screw_axes, joint_positions)) * Vs
        i += 1
        Tsb = FKinSpace(home_config, screw_axes, joint_positions)
        Vs = Adjoint(Tsb) * se3ToVec(MatrixLog6(TransInv(Tsb) * target_config))
        err = LA.norm(Vs[1:3]) > angular_tolerance || LA.norm(Vs[4:6]) > linear_tolerance
    end
    joint_positions, !err
end

# """
# *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
# """

"""
    ad(V)

Calculate the 6x6 matrix [adV] of the given 6-vector.

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
    ωmat = VecToso3(V[1:3])
    vcat(hcat(ωmat, zeros(3, 3)), hcat(VecToso3(V[4:6]), ωmat))
end

"""
    InverseDynamics(joint_positions, joint_velocities, joint_accelerations, gravity, tip_wrench, link_frames, spatial_inertias, screw_axes)

Computes inverse dynamics in the space frame for an open chain robot.
"""
function InverseDynamics(
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
    AdTi[n+1] = Adjoint(TransInv(link_frames[n+1]))
    Fi = copy(tip_wrench)
    joint_torques = zeros(eltype(joint_positions), n)

    for i = 1:n
        Mi *= link_frames[i]
        Ai[:, i] = Adjoint(TransInv(Mi)) * screw_axes[:, i]
        AdTi[i] = Adjoint(
            MatrixExp6(VecTose3(Ai[:, i] * -joint_positions[i])) * TransInv(link_frames[i]),
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
    MassMatrix(joint_positions, link_frames, spatial_inertias, screw_axes)

Computes the mass matrix of an open chain robot based on the given configuration.
"""
function MassMatrix(
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
        M[:, i] = InverseDynamics(
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
    VelQuadraticForces(joint_positions, joint_velocities, link_frames, spatial_inertias, screw_axes)

Computes the Coriolis and centripetal terms in the inverse dynamics of an open chain robot.
"""
function VelQuadraticForces(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    InverseDynamics(
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
    GravityForces(joint_positions, gravity, link_frames, spatial_inertias, screw_axes)

Computes the joint forces/torques an open chain robot requires to overcome gravity at its configuration.
"""
function GravityForces(
    joint_positions::AbstractVector,
    gravity::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    InverseDynamics(
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
    EndEffectorForces(joint_positions, tip_wrench, link_frames, spatial_inertias, screw_axes)

Computes the joint forces/torques an open chain robot requires only to create the end-effector force `tip_wrench`.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `tip_wrench`: the spatial force applied by the end-effector expressed in frame `{n+1}`.
- `link_frames`: the list of link frames `i` relative to `i-1` at the home position.
- `spatial_inertias`: the spatial inertia matrices `Gi` of the links.
- `screw_axes`: the screw axes `Si` of the joints in a space frame, in the format of a matrix with axes as the columns.

Returns the joint forces and torques required only to create the end-effector force `tip_wrench`.
This function calls InverseDynamics with `gravity = 0`, `joint_velocities = 0`, and `joint_accelerations = 0`.

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
 [0.010267 0.0 … 0.0 0.0; 0.0 0.010267 … 0.0 0.0; … ; 0.0 0.0 … 3.7 0.0; 0.0 0.0 … 0.0 3.7]
 [0.22689 0.0 … 0.0 0.0; 0.0 0.22689 … 0.0 0.0; … ; 0.0 0.0 … 8.393 0.0; 0.0 0.0 … 0.0 8.393]
 [0.0494433 0.0 … 0.0 0.0; 0.0 0.0494433 … 0.0 0.0; … ; 0.0 0.0 … 2.275 0.0; 0.0 0.0 … 0.0 2.275]

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

julia> EndEffectorForces(joint_positions, tip_wrench, link_frames, spatial_inertias, screw_axes)
3-element Vector{Float64}:
 1.4095460782639782
 1.857714972318063
 1.392409
```
"""
function EndEffectorForces(
    joint_positions::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    InverseDynamics(
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
    ForwardDynamics(joint_positions, joint_velocities, joint_torques, gravity, tip_wrench, link_frames, spatial_inertias, screw_axes)

Computes forward dynamics in the space frame for an open chain robot.
"""
function ForwardDynamics(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_torques::AbstractVector,
    gravity::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    LA.inv(MassMatrix(joint_positions, link_frames, spatial_inertias, screw_axes)) * (
        joint_torques - VelQuadraticForces(
            joint_positions,
            joint_velocities,
            link_frames,
            spatial_inertias,
            screw_axes,
        ) -
        GravityForces(joint_positions, gravity, link_frames, spatial_inertias, screw_axes) -
        EndEffectorForces(
            joint_positions,
            tip_wrench,
            link_frames,
            spatial_inertias,
            screw_axes,
        )
    )
end

"""
    EulerStep(joint_positions, joint_velocities, joint_accelerations, timestep)

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
julia> EulerStep([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], [2, 1.5, 1], 0.1)
([0.11000000000000001, 0.12000000000000001, 0.13], [0.30000000000000004, 0.35000000000000003, 0.4])
```
"""
function EulerStep(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_accelerations::AbstractVector,
    timestep::Number,
)
    joint_positions + timestep * joint_velocities,
    joint_velocities + timestep * joint_accelerations
end

"""
    InverseDynamicsTrajectory(joint_position_traj, joint_velocity_traj, joint_acceleration_traj, gravity, tip_wrench_traj, link_frames, spatial_inertias, screw_axes)

Calculates the joint forces/torques required to move the serial chain along the given trajectory using inverse dynamics.
"""
function InverseDynamicsTrajectory(
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

    for i = 1:size(joint_position_traj, 2)
        joint_torque_traj[:, i] = InverseDynamics(
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
    ForwardDynamicsTrajectory(joint_positions, joint_velocities, joint_torque_traj, gravity, tip_wrench_traj, link_frames, spatial_inertias, screw_axes, timestep, integration_resolution)

Simulates the motion of a serial chain given an open-loop history of joint forces/torques.
"""
function ForwardDynamicsTrajectory(
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

    for i = 1:size(joint_torque_traj, 2)-1
        for j = 1:integration_resolution
            joint_accelerations = ForwardDynamics(
                joint_positions,
                joint_velocities,
                joint_torque_traj[:, i],
                gravity,
                tip_wrench_traj[:, i],
                link_frames,
                spatial_inertias,
                screw_axes,
            )
            joint_positions, joint_velocities = EulerStep(
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
    CubicTimeScaling(total_time, t)

Computes s(t) for a cubic time scaling.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> CubicTimeScaling(2, 0.6)
0.21600000000000003
```
"""
CubicTimeScaling(total_time::Number, t::Number) = 3(t / total_time)^2 - 2(t / total_time)^3

"""
    QuinticTimeScaling(total_time, t)

Computes s(t) for a quintic time scaling.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> QuinticTimeScaling(2, 0.6)
0.16308
```
"""
QuinticTimeScaling(total_time::Number, t::Number) =
    10(t / total_time)^3 - 15(t / total_time)^4 + 6(t / total_time)^5

"""
    JointTrajectory(joint_position_start, joint_position_end, total_time, N, method)

Computes a straight-line trajectory in joint space.
"""
function JointTrajectory(
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
            s = CubicTimeScaling(total_time, timegap * (i - 1))
        else
            s = QuinticTimeScaling(total_time, timegap * (i - 1))
        end

        traj[:, i] = s * joint_position_end + (1 - s) * joint_position_start
    end

    traj'
end

"""
    ScrewTrajectory(transform_start, transform_end, total_time, N, method)

Computes a trajectory as a list of N SE(3) matrices corresponding to the screw motion about a space screw axis.
"""
function ScrewTrajectory(
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
            s = CubicTimeScaling(total_time, timegap * (i - 1))
        else
            s = QuinticTimeScaling(total_time, timegap * (i - 1))
        end

        traj[i] =
            transform_start *
            MatrixExp6(MatrixLog6(TransInv(transform_start) * transform_end) * s)
    end

    return traj
end

"""
    CartesianTrajectory(transform_start, transform_end, total_time, N, method)

Computes a trajectory as a list of N SE(3) matrices corresponding to the origin of the end-effector frame following a straight line.
"""
function CartesianTrajectory(
    transform_start::AbstractMatrix,
    transform_end::AbstractMatrix,
    total_time::Number,
    N::Integer,
    method::Integer,
)
    timegap = total_time / (N - 1)
    traj = Array{Array{Float64,2}}(undef, N)

    Rstart, pstart = TransToRp(transform_start)
    Rend, pend = TransToRp(transform_end)

    for i = 1:N
        if method == 3
            s = CubicTimeScaling(total_time, timegap * (i - 1))
        else
            s = QuinticTimeScaling(total_time, timegap * (i - 1))
        end

        traj[i] = vcat(
            hcat(
                Rstart * MatrixExp3(MatrixLog3(Rstart' * Rend) * s),
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
    ComputedTorque(joint_positions, joint_velocities, error_integral, gravity, link_frames, spatial_inertias, screw_axes, desired_joint_positions, desired_joint_velocities, desired_joint_accelerations, Kp, Ki, Kd)

Computes the joint control torques at a particular time instant.
"""
function ComputedTorque(
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
    MassMatrix(joint_positions, link_frames, spatial_inertias, screw_axes) * (
        Kp * e +
        Ki * (error_integral + e) +
        Kd * (desired_joint_velocities - joint_velocities)
    ) + InverseDynamics(
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
    SimulateControl(joint_positions, joint_velocities, gravity, tip_wrench_traj, link_frames, spatial_inertias,
                    screw_axes, desired_joint_position_traj, desired_joint_velocity_traj, desired_joint_acceleration_traj, estimated_gravity,
                    estimated_link_frames, estimated_spatial_inertias, Kp, Ki, Kd, timestep, integration_resolution)

Simulates the computed torque controller over a given desired trajectory.
"""
function SimulateControl(
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
        joint_torques = ComputedTorque(
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
            joint_accelerations = ForwardDynamics(
                current_positions,
                current_velocities,
                joint_torques,
                gravity,
                tip_wrench_traj[:, i],
                link_frames,
                spatial_inertias,
                screw_axes,
            )
            current_positions, current_velocities = EulerStep(
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

end # module

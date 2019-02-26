module ModernRoboticsBook

import LinearAlgebra
const linalg = LinearAlgebra

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
3-element Array{Float64,1}:
 0.2672612419124244
 0.5345224838248488
 0.8017837257372732
```
"""
Normalize(V::Array) = V / linalg.norm(V)

# """
# *** CHAPTER 3: RIGID-BODY MOTIONS ***
# """

"""
    RotInv(R)

Inverts a rotation matrix.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> RotInv([0 0 1; 1 0 0; 0 1 0])
3×3 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
 0  1  0
 0  0  1
 1  0  0
```
"""
RotInv(R::Array) = R'

"""
    VecToso3(omg)

Converts a 3-vector to an so(3) representation.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> VecToso3([1 2 3])
3×3 Array{Int64,2}:
  0  -3   2
  3   0  -1
 -2   1   0
```
"""
function VecToso3(omg::Array)
    [0      -omg[3]  omg[2];
     omg[3]       0 -omg[1];
     -omg[2] omg[1]       0;]
end

"""
    so3ToVec(so3mat)

Converts an so(3) representation to a 3-vector.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> so3ToVec([0 -3 2; 3 0 -1; -2 1 0])
3-element Array{Int64,1}:
 1
 2
 3
```
"""
function so3ToVec(so3mat::Array)
    [so3mat[3, 2], so3mat[1, 3], so3mat[2, 1]]
end

"""
    AxisAng3(expc3)

Converts a 3-vector of exponential coordinates for rotation into axis-angle form.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> AxisAng3([1, 2, 3])
([0.267261, 0.534522, 0.801784], 3.7416573867739413)
```
"""
AxisAng3(expc3::Array) = Normalize(expc3), linalg.norm(expc3)

"""
    MatrixExp3(so3mat)

Computes the matrix exponential of a matrix in so(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> MatrixExp3([0 -3 2; 3 0 -1; -2 1 0])
3×3 Array{Float64,2}:
 -0.694921   0.713521  0.0892929
 -0.192007  -0.303785  0.933192 
  0.692978   0.63135   0.348107 
```
"""
function MatrixExp3(so3mat::Array)
    omgtheta = so3ToVec(so3mat)
    if NearZero(linalg.norm(omgtheta))
        return linalg.I
    else
        θ = AxisAng3(omgtheta)[2]
        omgmat = so3mat / θ
        return linalg.I + sin(θ) * omgmat + (1 - cos(θ)) * omgmat * omgmat
    end
end

"""
    MatrixLog3(R)

Computes the matrix logarithm of a rotation matrix.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> MatrixLog3([0 0 1; 1 0 0; 0 1 0])
3×3 Array{Float64,2}:
  0.0     -1.2092   1.2092
  1.2092   0.0     -1.2092
 -1.2092   1.2092   0.0   
```
"""
function MatrixLog3(R::Array)
    acosinput = (linalg.tr(R) - 1) / 2.0
    if acosinput >= 1
        return zeros(3, 3)
    elseif acosinput <= -1
        if !NearZero(1 + R[3, 3])
            omg = (1.0 / √(2 * (1 + R[3, 3]))) * [R[1, 3], R[2, 3], 1 + R[3, 3]]
        elseif !NearZero(1 + R[2, 2])
            omg = (1.0 / √(2 * (1 + R[2, 2]))) * [R[1, 2], 1 + R[2, 2], R[3, 2]]
        else
            omg = (1.0 / √(2 * (1 + R[1, 1]))) * [1 + R[1, 1], R[2, 1], R[3, 1]]
        end
        return VecToso3(π * omg)
    else
        θ = acos(acosinput)
        return θ / 2.0 / sin(θ) * (R - R')
    end
end

"""
    RpToTrans(R, p)

Converts a rotation matrix and a position vector into homogeneous transformation matrix.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> RpToTrans([1 0 0; 0 0 -1; 0 1 0], [1, 2, 5])
4×4 Array{Int64,2}:
 1  0   0  1
 0  0  -1  2
 0  1   0  5
 0  0   0  1
```
"""
RpToTrans(R::Array, p::Array) = vcat(hcat(R, p), [0 0 0 1])

"""
    TransToRp(T)

Converts a homogeneous transformation matrix into a rotation matrix and position vector.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> TransToRp([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
([1 0 0; 0 0 -1; 0 1 0], [0, 0, 3])
```
"""
TransToRp(T::Array) = T[1:3, 1:3], T[1:3, 4]

"""
    TransInv(T)

Inverts a homogeneous transformation matrix.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> TransInv([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
4×4 Array{Int64,2}:
 1   0  0   0
 0   0  1  -3
 0  -1  0   0
 0   0  0   1
```
"""
function TransInv(T::Array)
    R, p = TransToRp(T)
    vcat(hcat(R', -R' * p), [0 0 0 1])
end

"""
    VecTose3(V)

Converts a spatial velocity vector into a 4x4 matrix in se3.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> VecTose3([1 2 3 4 5 6])
4×4 Array{Float64,2}:
  0.0  -3.0   2.0  4.0
  3.0   0.0  -1.0  5.0
 -2.0   1.0   0.0  6.0
  0.0   0.0   0.0  0.0
```
"""
VecTose3(V::Array) = vcat(hcat(VecToso3(V[1:3]), V[4:6]), zeros(1, 4))

"""
    se3ToVec(se3mat)

Converts an se3 matrix into a spatial velocity vector.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> se3ToVec([0 -3 2 4; 3 0 -1 5; -2 1 0 6; 0 0 0 0])
6-element Array{Int64,1}:
 1
 2
 3
 4
 5
 6
```
"""
se3ToVec(se3mat::Array) = vcat([se3mat[3, 2],
                                se3mat[1, 3],
                                se3mat[2, 1]], se3mat[1:3, 4])

"""
    Adjoint(T)

Computes the adjoint representation of a homogeneous transformation matrix.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> Adjoint([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
6×6 Array{Float64,2}:
 1.0  0.0   0.0  0.0  0.0   0.0
 0.0  0.0  -1.0  0.0  0.0   0.0
 0.0  1.0   0.0  0.0  0.0   0.0
 0.0  0.0   3.0  1.0  0.0   0.0
 3.0  0.0   0.0  0.0  0.0  -1.0
 0.0  0.0   0.0  0.0  1.0   0.0
```
"""
function Adjoint(T::Array)
    R, p = TransToRp(T)
    vcat(hcat(R, zeros(3, 3)), hcat(VecToso3(p) * R, R))
end

"""
    ScrewToAxis(q, s, h)

Takes a parametric description of a screw axis and converts it to a normalized screw axis.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> ScrewToAxis([3; 0; 0], [0; 0; 1], 2)
6-element Array{Int64,1}:
  0
  0
  1
  0
 -3
  2
```
"""
ScrewToAxis(q::Array, s::Array, h::Number) = vcat(s, linalg.cross(q, s) + h * s)

"""
    AxisAng6(expc6)

Converts a 6-vector of exponential coordinates into screw axis-angle form.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> AxisAng6([1, 0, 0, 1, 2, 3])
([1.0, 0.0, 0.0, 1.0, 2.0, 3.0], 1.0)
```
"""
function AxisAng6(expc6::Array)
    θ = linalg.norm(expc6[1:3])
    if NearZero(θ)
        θ = linalg.norm(expc6[3:6])
    end
    expc6 / θ, θ
end

"""
    MatrixExp6(se3mat)

Computes the matrix exponential of an se3 representation of exponential coordinates.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> MatrixExp6([0 0 0 0; 0 0 -1.57079632 2.35619449; 0 1.57079632 0 2.35619449; 0 0 0 0])
4×4 Array{Float64,2}:
 1.0  0.0         0.0        0.0       
 0.0  6.7949e-9  -1.0        1.01923e-8
 0.0  1.0         6.7949e-9  3.0       
 0.0  0.0         0.0        1.0       
```
"""
function MatrixExp6(se3mat::Array)
    omgtheta = so3ToVec(se3mat[1:3, 1:3])
    if NearZero(linalg.norm(omgtheta))
        return vcat(hcat(linalg.I, se3mat[1:3, 4]), [0 0 0 1])
    else
        θ = AxisAng3(omgtheta)[2]
        omgmat = se3mat[1:3, 1:3] / θ
        return vcat(hcat(MatrixExp3(se3mat[1:3, 1:3]),
                         (linalg.I * θ +
                          (1 - cos(θ)) * omgmat +
                          (θ - sin(θ)) * omgmat * omgmat) *
                         se3mat[1:3, 4] / θ),
                    [0 0 0 1])
    end
end

"""
    MatrixLog6(T)

Computes the matrix logarithm of a homogeneous transformation matrix.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> MatrixLog6([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
4×4 Array{Float64,2}:
 0.0  0.0      0.0     0.0    
 0.0  0.0     -1.5708  2.35619
 0.0  1.5708   0.0     2.35619
 0.0  0.0      0.0     0.0    
```
"""
function MatrixLog6(T::Array)
    R, p = TransToRp(T)
    omgmat = MatrixLog3(R)
    if omgmat == zeros(3, 3)
        return vcat(hcat(zeros(3, 3), T[1:3, 4]), [0 0 0 0])
    else
        θ = acos((linalg.tr(R) - 1) / 2.0)
        return vcat(hcat(omgmat,
                         (linalg.I - omgmat / 2.0 +
                          (1.0 / θ - 1.0 / tan(θ / 2.0) / 2) *
                          omgmat * omgmat / θ) * T[1:3, 4]),
                    [0 0 0 0])
    end
end

"""
    ProjectToSO3(mat)

Returns a projection of mat into SO(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> ProjectToSO3([0.675 0.150  0.720; 0.370 0.771 -0.511; -0.630 0.619  0.472])
3×3 Array{Float64,2}:
  0.679011  0.148945   0.718859
  0.373207  0.773196  -0.512723
 -0.632187  0.616428   0.469421
```
"""
function ProjectToSO3(mat::Array)
    F  = linalg.svd(mat)
    R = F.U * F.Vt
    if linalg.det(R) < 0
        # In this case the result may be far from mat.
        R[:, F.s[3, 3]] = -R[:, F.s[3, 3]]
    end
    return R
end

"""
    ProjectToSE3(mat)

Returns a projection of mat into SE(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> ProjectToSE3([0.675 0.150 0.720 1.2; 0.370 0.771 -0.511 5.4; -0.630 0.619 0.472 3.6; 0.003 0.002 0.010 0.9])
4×4 Array{Float64,2}:
  0.679011  0.148945   0.718859  1.2
  0.373207  0.773196  -0.512723  5.4
 -0.632187  0.616428   0.469421  3.6
  0.0       0.0        0.0       1.0
```
"""
ProjectToSE3(mat::Array) = RpToTrans(ProjectToSO3(mat[1:3, 1:3]), mat[1:3, 4])

"""
    DistanceToSO3(mat)

Returns the Frobenius norm to describe the distance of mat from the SO(3) manifold.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> DistanceToSO3([1.0 0.0 0.0; 0.0 0.1 -0.95; 0.0 1.0 0.1])
0.08835298523536149
```
"""
DistanceToSO3(mat::Array) = linalg.det(mat) > 0 ? linalg.norm(mat'mat - linalg.I) : 1e+9

"""
    DistanceToSE3(mat)

Returns the Frobenius norm to describe the distance of mat from the SE(3) manifold.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> DistanceToSE3([1.0 0.0 0.0 1.2; 0.0 0.1 -0.95 1.5; 0.0 1.0 0.1 -0.9; 0.0 0.0 0.1 0.98])
0.13493053768513638
```
"""
function DistanceToSE3(mat::Array)
    matR = mat[1:3, 1:3]
    if linalg.det(matR) > 0
        linalg.norm(hcat(vcat(matR'matR, zeros(1, 3)), mat[4, :]) - linalg.I)
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
TestIfSO3(mat::Array) = abs(DistanceToSO3(mat)) < 1e-3

"""
    TestIfSE3(mat)

Returns true if mat is close to or on the manifold SE(3).

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> TestIfSE3([1.0 0.0 0.0 1.2; 0.0 0.1 -0.95 1.5; 0.0 1.0 0.1 -0.9; 0.0 0.0 0.1 0.98])
false
```
"""
TestIfSE3(mat::Array) = abs(DistanceToSE3(mat)) < 1e-3

# """
# *** CHAPTER 4: FORWARD KINEMATICS ***
# """

"""
    FKinBody(M, Blist, thetalist)

Computes forward kinematics in the body frame for an open chain robot.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> M = [ -1  0  0  0 ;
              0  1  0  6 ;
              0  0 -1  2 ;
              0  0  0  1 ];

julia> Blist = [  0  0 -1  2  0  0   ;
                  0  0  0  0  1  0   ;
                  0  0  1  0  0  0.1 ]';

julia> thetalist = [ π/2, 3, π ];

julia> FKinBody(M, Blist, thetalist)
4×4 Array{Float64,2}:
 -1.14424e-17  1.0           0.0  -5.0    
  1.0          1.14424e-17   0.0   4.0    
  0.0          0.0          -1.0   1.68584
  0.0          0.0           0.0   1.0    
```
"""
function FKinBody(M::AbstractMatrix, Blist::AbstractMatrix, thetalist::Array)
    for i = 1:length(thetalist)
        M *= MatrixExp6(VecTose3(Blist[:, i] * thetalist[i]))
    end
    M
end

"""
    FKinSpace(M, Slist, thetalist)

Computes forward kinematics in the space frame for an open chain robot.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> M = [ -1  0  0  0 ;
              0  1  0  6 ;
              0  0 -1  2 ;
              0  0  0  1 ];

julia> Slist = [  0  0  1  4  0  0   ;
                  0  0  0  0  1  0   ;
                  0  0 -1 -6  0 -0.1 ]';

julia> thetalist = [ π/2, 3, π ];

julia> FKinSpace(M, Slist, thetalist)
4×4 Array{Float64,2}:
 -1.14424e-17  1.0           0.0  -5.0    
  1.0          1.14424e-17   0.0   4.0    
  0.0          0.0          -1.0   1.68584
  0.0          0.0           0.0   1.0    
```
"""
function FKinSpace(M::AbstractMatrix, Slist::AbstractMatrix, thetalist::Array)
    for i = length(thetalist):-1:1
        M = MatrixExp6(VecTose3(Slist[:, i] * thetalist[i])) * M
    end
    M
end

# """
# *** CHAPTER 5: VELOCITY KINEMATICS AND STATICS ***
# """

"""
    JacobianBody(Blist, thetalist)

Computes the body Jacobian for an open chain robot.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> Blist = [0 0 1   0 0.2 0.2;
                1 0 0   2   0   3;
                0 1 0   0   2   1;
                1 0 0 0.2 0.3 0.4]';

julia> thetalist = [0.2, 1.1, 0.1, 1.2];

julia> JacobianBody(Blist, thetalist)
6×4 Array{Float64,2}:
 -0.0452841  0.995004    0.0       1.0
  0.743593   0.0930486   0.362358  0.0
 -0.667097   0.0361754  -0.932039  0.0
  2.32586    1.66809     0.564108  0.2
 -1.44321    2.94561     1.43307   0.3
 -2.0664     1.82882    -1.58869   0.4
```
"""
function JacobianBody(Blist::AbstractMatrix, thetalist::Array)
    T = linalg.I
    Jb = copy(Blist)
    for i = length(thetalist)-1:-1:1
        T *= MatrixExp6(VecTose3(Blist[:, i+1] * -thetalist[i+1]))
        Jb[:, i] = Adjoint(T) * Blist[:, i]
    end
    Jb
end

"""
    JacobianSpace(Slist, thetalist)

Computes the space Jacobian for an open chain robot.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> Slist = [0 0 1   0 0.2 0.2;
                1 0 0   2   0   3;
                0 1 0   0   2   1;
                1 0 0 0.2 0.3 0.4]';

julia> thetalist = [0.2, 1.1, 0.1, 1.2];

julia> JacobianSpace(Slist, thetalist)
6×4 Array{Float64,2}:
 0.0  0.980067  -0.0901156   0.957494 
 0.0  0.198669   0.444554    0.284876 
 1.0  0.0        0.891207   -0.0452841
 0.0  1.95219   -2.21635    -0.511615 
 0.2  0.436541  -2.43713     2.77536  
 0.2  2.96027    3.23573     2.22512  
```
"""
function JacobianSpace(Slist::AbstractMatrix, thetalist::Array)
    T = linalg.I
    Js = copy(Slist)
    for i = 2:length(thetalist)
        T *= MatrixExp6(VecTose3(Slist[:, i - 1] * thetalist[i - 1]))
        Js[:, i] = Adjoint(T) * Slist[:, i]
    end
    Js
end

# """
# *** CHAPTER 6: INVERSE KINEMATICS ***
# """

"""
    IKinBody(Blist, M, T, thetalist0, eomg, ev)

Computes inverse kinematics in the body frame for an open chain robot.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> Blist = [  0  0 -1  2  0  0   ;
                  0  0  0  0  1  0   ;
                  0  0  1  0  0  0.1 ]';

julia> M = [ -1  0  0  0 ;
              0  1  0  6 ;
              0  0 -1  2 ;
              0  0  0  1 ];

julia> T = [ 0  1  0     -5 ;
             1  0  0      4 ;
             0  0 -1 1.6858 ;
             0  0  0      1 ];

julia> thetalist0 = [1.5, 2.5, 3];

julia> eomg, ev = 0.01, 0.001;

julia> IKinBody(Blist, M, T, thetalist0, eomg, ev)
([1.57074, 2.99967, 3.14154], true)
```
"""
function IKinBody(Blist::AbstractMatrix,
                      M::AbstractMatrix,
                      T::AbstractMatrix,
             thetalist0::Array,
                   eomg::Number,
                     ev::Number)
    thetalist = copy(thetalist0)
    i = 0
    maxiterations = 20
    Vb = se3ToVec(MatrixLog6(TransInv(FKinBody(M, Blist, thetalist)) * T))
    err = linalg.norm(Vb[1:3]) > eomg || linalg.norm(Vb[4:6]) > ev
    while err && i < maxiterations
        thetalist += linalg.pinv(JacobianBody(Blist, thetalist)) * Vb
        i += 1
        Vb = se3ToVec(MatrixLog6(TransInv(FKinBody(M, Blist, thetalist)) * T))
        err = linalg.norm(Vb[1:3]) > eomg || linalg.norm(Vb[4:6]) > ev
    end
    return thetalist, !err
end

"""
    IKinSpace(Slist, M, T, thetalist0, eomg, ev)

Computes inverse kinematics in the space frame for an open chain robot.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> Slist = [  0  0  1  4  0  0   ;
                  0  0  0  0  1  0   ;
                  0  0 -1 -6  0 -0.1 ]';

julia> M = [ -1  0  0  0 ;
              0  1  0  6 ;
              0  0 -1  2 ;
              0  0  0  1 ];

julia> T = [ 0  1  0     -5 ;
             1  0  0      4 ;
             0  0 -1 1.6858 ;
             0  0  0      1 ];

julia> thetalist0 = [1.5, 2.5, 3];

julia> eomg, ev = 0.01, 0.001;

julia> IKinSpace(Slist, M, T, thetalist0, eomg, ev)
([1.57074, 2.99966, 3.14153], true)
```
"""
function IKinSpace(Slist::AbstractMatrix,
                       M::AbstractMatrix,
                       T::AbstractMatrix,
              thetalist0::Array,
                    eomg::Number,
                      ev::Number)
    thetalist = copy(thetalist0)
    i = 0
    maxiterations = 20
    Tsb = FKinSpace(M, Slist, thetalist)
    Vs = Adjoint(Tsb) * se3ToVec(MatrixLog6(TransInv(Tsb) * T))
    err = linalg.norm(Vs[1:3]) > eomg || linalg.norm(Vs[4:6]) > ev
    while err && i < maxiterations
        thetalist += linalg.pinv(JacobianSpace(Slist, thetalist)) * Vs
        i += 1
        Tsb = FKinSpace(M, Slist, thetalist)
        Vs = Adjoint(Tsb) * se3ToVec(MatrixLog6(TransInv(Tsb) * T))
        err = linalg.norm(Vs[1:3]) > eomg || linalg.norm(Vs[4:6]) > ev
    end
    thetalist, !err
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
6×6 Array{Float64,2}:
  0.0  -3.0   2.0   0.0   0.0   0.0
  3.0   0.0  -1.0   0.0   0.0   0.0
 -2.0   1.0   0.0   0.0   0.0   0.0
  0.0  -6.0   5.0   0.0  -3.0   2.0
  6.0   0.0  -4.0   3.0   0.0  -1.0
 -5.0   4.0   0.0  -2.0   1.0   0.0
```
"""
function ad(V::Array)
    omgmat = VecToso3(V[1:3])
    vcat(hcat(omgmat, zeros(3, 3)),
         hcat(VecToso3(V[4:6]), omgmat))
end

"""
    InverseDynamics(thetalist, dthetalist, ddthetalist, g, Ftip, Mlist, Glist, Slist)

Computes inverse dynamics in the space frame for an open chain robot.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> InverseDynamics(thetalist, dthetalist, ddthetalist, g, Ftip, Mlist, Glist, Slist)
3-element Array{Float64,1}:
  74.69616155287451 
 -33.06766015851458 
  -3.230573137901424
```
"""
function InverseDynamics(thetalist::Array,
                        dthetalist::Array,
                       ddthetalist::Array,
                                 g::Array,
                              Ftip::Array,
                             Mlist::Array,
                             Glist::Array,
                             Slist::AbstractMatrix)
    n = length(thetalist)
    Mi = linalg.I
    Ai = zeros(6, n)
    AdTi = Array{Array{Float64, 2}}(undef, n + 1)
    Vi = zeros(6, n + 1)
    Vdi = zeros(6, n + 1)
    Vdi[:, 1] = vcat([0, 0, 0], -g)
    AdTi[n+1] = Adjoint(TransInv(Mlist[n+1]))
    Fi = copy(Ftip)
    taulist = zeros(n)

    for i = 1:n
        Mi *= Mlist[i]
        Ai[:, i] = Adjoint(TransInv(Mi)) * Slist[:, i]
        AdTi[i] = Adjoint(MatrixExp6(VecTose3(Ai[:, i] * -thetalist[i])) *
                          TransInv(Mlist[i]))
        Vi[:, i + 1] = AdTi[i] * Vi[:,i] + Ai[:, i] * dthetalist[i]
        Vdi[:, i + 1] = AdTi[i] * Vdi[:, i] + Ai[:, i] * ddthetalist[i] +
                        ad(Vi[:, i + 1]) * Ai[:, i] * dthetalist[i]
    end

    for i = n:-1:1
        Fi = AdTi[i + 1]' * Fi + Glist[i] * Vdi[:, i + 1] -
             ad(Vi[:, i + 1])' * Glist[i] * Vi[:, i + 1]
        taulist[i] = Fi' * Ai[:, i]
    end

    return taulist
end

"""
    MassMatrix(thetalist, Mlist, Glist, Slist)

Computes the mass matrix of an open chain robot based on the given configuration.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> MassMatrix(thetalist, Mlist, Glist, Slist)
3×3 Array{Float64,2}:
 22.5433      -0.307147  -0.00718426
 -0.307147     1.96851    0.432157  
 -0.00718426   0.432157   0.191631  
```
"""
function MassMatrix(thetalist::Array,
                        Mlist::Array,
                        Glist::Array,
                        Slist::AbstractMatrix)
    n = length(thetalist)
    M = zeros(n, n)

    for i = 1:n
        ddthetalist = zeros(n)
        ddthetalist[i] = 1
        M[:, i] = InverseDynamics(thetalist, zeros(n), ddthetalist, [0, 0, 0],
                                  [0, 0, 0, 0, 0, 0], Mlist, Glist, Slist)
    end

    return M
end

"""
    VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist)

Computes the Coriolis and centripetal terms in the inverse dynamics of an open chain robot.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist)
3-element Array{Float64,1}:
  0.26453118054501235 
 -0.0550515682891655  
 -0.006891320068248911
```
"""
function VelQuadraticForces(thetalist::Array,
                           dthetalist::Array,
                                Mlist::Array,
                                Glist::Array,
                                Slist::AbstractMatrix)
    InverseDynamics(thetalist, dthetalist, zeros(length(thetalist)),
                    [0, 0, 0], [0, 0, 0, 0, 0, 0], Mlist, Glist, Slist)
end

"""
    GravityForces(thetalist, g, Mlist, Glist, Slist)

Computes the joint forces/torques an open chain robot requires to overcome gravity at its configuration.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> GravityForces(thetalist, g, Mlist, Glist, Slist)
3-element Array{Float64,1}:
  28.40331261821983  
 -37.64094817177068  
  -5.4415891999683605
```
"""
function GravityForces(thetalist::Array,
                               g::Array,
                           Mlist::Array,
                           Glist::Array,
                           Slist::AbstractMatrix)
    n = length(thetalist)
    InverseDynamics(thetalist, zeros(n), zeros(n), g, [0, 0, 0, 0, 0, 0], Mlist, Glist, Slist)
end

"""
    EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist)

Computes the joint forces/torques an open chain robot requires only to create the end-effector force `Ftip`.

# Arguments
- `thetalist`: the ``n``-vector of joint variables.
- `Ftip`: the spatial force applied by the end-effector expressed in frame `{n+1}`.
- `Mlist`: the list of link frames `i` relative to `i-1` at the home position.
- `Glist`: the spatial inertia matrices `Gi` of the links.
- `Slist`: the screw axes `Si` of the joints in a space frame, in the format of a matrix with axes as the columns.

Returns the joint forces and torques required only to create the end-effector force `Ftip`.
This function calls InverseDynamics with `g = 0`, `dthetalist = 0`, and `ddthetalist = 0`.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> import LinearAlgebra

julia> const linalg = LinearAlgebra;

julia> thetalist = [0.1, 0.1, 0.1]
3-element Array{Float64,1}:
 0.1
 0.1
 0.1

julia> Ftip = [1, 1, 1, 1, 1, 1]
6-element Array{Int64,1}:
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
4×4 Array{Float64,2}:
 1.0  0.0  0.0  0.0     
 0.0  1.0  0.0  0.0     
 0.0  0.0  1.0  0.089159
 0.0  0.0  0.0  1.0     

julia> M12 = [ 0 0 1    0.28;
               0 1 0 0.13585;
              -1 0 0       0;
               0 0 0       1]
4×4 Array{Float64,2}:
  0.0  0.0  1.0  0.28   
  0.0  1.0  0.0  0.13585
 -1.0  0.0  0.0  0.0    
  0.0  0.0  0.0  1.0    

julia> M23 = [1 0 0       0;
              0 1 0 -0.1197;
              0 0 1   0.395;
              0 0 0       1]
4×4 Array{Float64,2}:
 1.0  0.0  0.0   0.0   
 0.0  1.0  0.0  -0.1197
 0.0  0.0  1.0   0.395 
 0.0  0.0  0.0   1.0   

julia> M34 = [1 0 0       0;
              0 1 0       0;
              0 0 1 0.14225;
              0 0 0       1]
4×4 Array{Float64,2}:
 1.0  0.0  0.0  0.0    
 0.0  1.0  0.0  0.0    
 0.0  0.0  1.0  0.14225
 0.0  0.0  0.0  1.0    

julia> Mlist = [M01, M12, M23, M34]
4-element Array{Array{Float64,2},1}:
 [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.089159; 0.0 0.0 0.0 1.0] 
 [0.0 0.0 1.0 0.28; 0.0 1.0 0.0 0.13585; -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0]
 [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 -0.1197; 0.0 0.0 1.0 0.395; 0.0 0.0 0.0 1.0]
 [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.14225; 0.0 0.0 0.0 1.0]  

julia> G1 = linalg.Diagonal([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
6×6 LinearAlgebra.Diagonal{Float64,Array{Float64,1}}:
 0.010267   ⋅         ⋅        ⋅    ⋅    ⋅ 
  ⋅        0.010267   ⋅        ⋅    ⋅    ⋅ 
  ⋅         ⋅        0.00666   ⋅    ⋅    ⋅ 
  ⋅         ⋅         ⋅       3.7   ⋅    ⋅ 
  ⋅         ⋅         ⋅        ⋅   3.7   ⋅ 
  ⋅         ⋅         ⋅        ⋅    ⋅   3.7

julia> G2 = linalg.Diagonal([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
6×6 LinearAlgebra.Diagonal{Float64,Array{Float64,1}}:
 0.22689   ⋅        ⋅          ⋅      ⋅      ⋅   
  ⋅       0.22689   ⋅          ⋅      ⋅      ⋅   
  ⋅        ⋅       0.0151074   ⋅      ⋅      ⋅   
  ⋅        ⋅        ⋅         8.393   ⋅      ⋅   
  ⋅        ⋅        ⋅          ⋅     8.393   ⋅   
  ⋅        ⋅        ⋅          ⋅      ⋅     8.393

julia> G3 = linalg.Diagonal([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])
6×6 LinearAlgebra.Diagonal{Float64,Array{Float64,1}}:
 0.0494433   ⋅          ⋅         ⋅      ⋅      ⋅   
  ⋅         0.0494433   ⋅         ⋅      ⋅      ⋅   
  ⋅          ⋅         0.004095   ⋅      ⋅      ⋅   
  ⋅          ⋅          ⋅        2.275   ⋅      ⋅   
  ⋅          ⋅          ⋅         ⋅     2.275   ⋅   
  ⋅          ⋅          ⋅         ⋅      ⋅     2.275

julia> Glist = [G1, G2, G3]
3-element Array{LinearAlgebra.Diagonal{Float64,Array{Float64,1}},1}:
 [0.010267 0.0 … 0.0 0.0; 0.0 0.010267 … 0.0 0.0; … ; 0.0 0.0 … 3.7 0.0; 0.0 0.0 … 0.0 3.7]      
 [0.22689 0.0 … 0.0 0.0; 0.0 0.22689 … 0.0 0.0; … ; 0.0 0.0 … 8.393 0.0; 0.0 0.0 … 0.0 8.393]    
 [0.0494433 0.0 … 0.0 0.0; 0.0 0.0494433 … 0.0 0.0; … ; 0.0 0.0 … 2.275 0.0; 0.0 0.0 … 0.0 2.275]

julia> Slist = [ 1  0  1      0  1      0;
                 0  1  0 -0.089  0      0;
                 0  1  0 -0.089  0  0.425]'
6×3 LinearAlgebra.Adjoint{Float64,Array{Float64,2}}:
 1.0   0.0     0.0  
 0.0   1.0     1.0  
 1.0   0.0     0.0  
 0.0  -0.089  -0.089
 1.0   0.0     0.0  
 0.0   0.0     0.425

julia> EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist)
3-element Array{Float64,1}:
 1.4095460782639782
 1.8577149723180628
 1.392409          
```
"""
function EndEffectorForces(thetalist::Array,
                                Ftip::Array,
                               Mlist::Array,
                               Glist::Array,
                               Slist::AbstractMatrix)
    n = length(thetalist)
    InverseDynamics(thetalist, zeros(n), zeros(n), [0, 0, 0], Ftip, Mlist, Glist, Slist)
end

"""
    ForwardDynamics(thetalist, dthetalist, taulist, g, Ftip, Mlist, Glist, Slist)

Computes forward dynamics in the space frame for an open chain robot.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> ForwardDynamics(thetalist, dthetalist, taulist, g, Ftip, Mlist, Glist, Slist)
3-element Array{Float64,1}:
  -0.9739290670855626
  25.584667840340558 
 -32.91499212478149  
```
"""
function ForwardDynamics(thetalist::Array,
                        dthetalist::Array,
                           taulist::Array,
                                 g::Array,
                              Ftip::Array,
                             Mlist::Array,
                             Glist::Array,
                             Slist::AbstractMatrix)
    linalg.inv(MassMatrix(thetalist, Mlist, Glist, Slist)) *
    (taulist - VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist)
             - GravityForces(thetalist, g, Mlist, Glist, Slist)
             - EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist))
end

"""
    EulerStep(thetalist, dthetalist, ddthetalist, dt)

Compute the joint angles and velocities at the next timestep using first order Euler integration.

# Arguments
- `thetalist`: the ``n``-vector of joint variables.
- `dthetalist`: the ``n``-vector of joint rates.
- `ddthetalist`: the ``n``-vector of joint accelerations.
- `dt`: the timestep delta t.

# Return
- `thetalistNext`: the vector of joint variables after `dt` from first order Euler integration.
- `dthetalistNext`: the vector of joint rates after `dt` from first order Euler integration.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> EulerStep([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], [2, 1.5, 1], 0.1)
([0.11, 0.12, 0.13], [0.3, 0.35, 0.4])
```
"""
function EulerStep(thetalist::Array, dthetalist::Array, ddthetalist::Array, dt::Number)
    thetalist + dt * dthetalist, dthetalist + dt * ddthetalist
end

"""
    InverseDynamicsTrajectory(thetamat, dthetamat, ddthetamat, g, Ftipmat, Mlist, Glist, Slist)

Calculates the joint forces/torques required to move the serial chain along the given trajectory using inverse dynamics.
"""
function InverseDynamicsTrajectory(thetamat::Array,
                                  dthetamat::Array,
                                 ddthetamat::Array,
                                          g::Array,
                                    Ftipmat::Array,
                                      Mlist::Array,
                                      Glist::Array,
                                      Slist::AbstractMatrix)
    thetamat = thetamat'
    dthetamat = dthetamat'
    ddthetamat = ddthetamat'
    Ftipmat = Ftipmat'
    taumat = copy(thetamat)

    for i = 1:size(thetamat, 2)
        taumat[:, i] = InverseDynamics(thetamat[:, i], dthetamat[:, i],
                                       ddthetamat[:, i], g, Ftipmat[:, i],
                                       Mlist, Glist, Slist)
    end

    taumat'
end

"""
    ForwardDynamicsTrajectory(thetalist, dthetalist, taumat, g, Ftipmat, Mlist, Glist, Slist, dt, intRes)

Simulates the motion of a serial chain given an open-loop history of joint forces/torques.
"""
function ForwardDynamicsTrajectory(thetalist::Array,
                                  dthetalist::Array,
                                      taumat::AbstractMatrix,
                                           g::Array,
                                     Ftipmat::Array,
                                       Mlist::Array,
                                       Glist::Array,
                                       Slist::AbstractMatrix,
                                          dt::Number,
                                      intRes::Number)
    taumat = taumat'
    Ftipmat = Ftipmat'
    thetamat = copy(taumat)
    thetamat[:, 1] = thetalist
    dthetamat = copy(taumat)
    dthetamat[:, 1] = dthetalist

    for i = 1:size(taumat, 2)-1
        for j = 1:intRes
            ddthetalist = ForwardDynamics(thetalist, dthetalist, taumat[:, i], g, Ftipmat[:, i], Mlist, Glist, Slist)
            thetalist, dthetalist = EulerStep(thetalist, dthetalist, ddthetalist, 1.0 * dt / intRes)
        end

        thetamat[:, i + 1] = thetalist
        dthetamat[:, i + 1] = dthetalist
    end

    thetamat = thetamat'
    dthetamat = dthetamat'
    return thetamat, dthetamat
end

# """
# *** CHAPTER 9: TRAJECTORY GENERATION ***
# """

"""
    CubicTimeScaling(Tf, t)

Computes s(t) for a cubic time scaling.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> CubicTimeScaling(2, 0.6)
0.21600000000000003
```
"""
CubicTimeScaling(Tf::Number, t::Number) = 3(t / Tf)^2 - 2(t / Tf)^3

"""
    QuinticTimeScaling(Tf, t)

Computes s(t) for a quintic time scaling.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> QuinticTimeScaling(2, 0.6)
0.16308
```
"""
QuinticTimeScaling(Tf::Number, t::Number) = 10(t / Tf)^3 - 15(t / Tf)^4 + 6(t / Tf)^5

"""
    JointTrajectory(thetastart, thetaend, Tf, N, method)

Computes a straight-line trajectory in joint space.
"""
function JointTrajectory(thetastart::Array, thetaend::Array, Tf::Number, N::Integer, method::Integer)
    timegap = Tf / (N - 1)
    traj = zeros(length(thetastart), N)

    for i = 1:N
        if method == 3
            s = CubicTimeScaling(Tf, timegap * (i - 1))
        else
            s = QuinticTimeScaling(Tf, timegap * (i - 1))
        end

        traj[:, i] = s * thetaend + (1 - s) * thetastart
    end

    traj'
end

"""
    ScrewTrajectory(Xstart, Xend, Tf, N, method)

Computes a trajectory as a list of N SE(3) matrices corresponding to the screw motion about a space screw axis.
"""
function ScrewTrajectory(Xstart::Array, Xend::Array, Tf::Number, N::Integer, method::Integer)
    timegap = Tf / (N - 1)
    traj = Array{Array{Float64, 2}}(undef, N)

    for i = 1:N
        if method == 3
            s = CubicTimeScaling(Tf, timegap * (i - 1))
        else
            s = QuinticTimeScaling(Tf, timegap * (i - 1))
        end

        traj[i] = Xstart * MatrixExp6(MatrixLog6(TransInv(Xstart) * Xend) * s)
    end

    return traj
end

"""
    CartesianTrajectory(Xstart, Xend, Tf, N, method)

Computes a trajectory as a list of N SE(3) matrices corresponding to the origin of the end-effector frame following a straight line.
"""
function CartesianTrajectory(Xstart::Array, Xend::Array, Tf::Number, N::Integer, method::Integer)
    timegap = Tf / (N - 1)
    traj = Array{Array{Float64, 2}}(undef, N)

    Rstart, pstart = TransToRp(Xstart)
    Rend, pend = TransToRp(Xend)

    for i = 1:N
        if method == 3
            s = CubicTimeScaling(Tf, timegap * (i - 1))
        else
            s = QuinticTimeScaling(Tf, timegap * (i - 1))
        end

        traj[i] = vcat(hcat(Rstart * MatrixExp3(MatrixLog3(Rstart' * Rend) * s), s * pend + (1 - s) * pstart), [0 0 0 1])
    end

    return traj
end

# """
# *** CHAPTER 11: ROBOT CONTROL ***
# """

"""
    ComputedTorque(thetalist, dthetalist, eint, g, Mlist, Glist, Slist, thetalistd, dthetalistd, ddthetalistd, Kp, Ki, Kd)

Computes the joint control torques at a particular time instant.
"""
function ComputedTorque(thetalist::Array,
                       dthetalist::Array,
                             eint::Array,
                                g::Array,
                            Mlist::Array,
                            Glist::Array,
                            Slist::AbstractMatrix,
                       thetalistd::Array,
                      dthetalistd::Array,
                     ddthetalistd::Array,
                               Kp::Number,
                               Ki::Number,
                               Kd::Number)
    e = thetalistd - thetalist
    MassMatrix(thetalist, Mlist, Glist, Slist) *
    (Kp * e + Ki * (eint + e) + Kd * (dthetalistd - dthetalist)) +
    InverseDynamics(thetalist, dthetalist, ddthetalistd, g,
                    [0, 0, 0, 0, 0, 0], Mlist, Glist, Slist)
end

"""
    SimulateControl(thetalist, dthetalist, g, Ftipmat, Mlist, Glist,
                    Slist, thetamatd, dthetamatd, ddthetamatd, gtilde,
                    Mtildelist, Gtildelist, Kp, Ki, Kd, dt, intRes)

Simulates the computed torque controller over a given desired trajectory.
"""
function SimulateControl(thetalist::Array,
                        dthetalist::Array,
                                 g::Array,
                           Ftipmat::Array,
                             Mlist::Array,
                             Glist::Array,
                             Slist::AbstractMatrix,
                         thetamatd::Array,
                        dthetamatd::Array,
                       ddthetamatd::Array,
                            gtilde::Array,
                        Mtildelist::Array,
                        Gtildelist::Array,
                                Kp::Number,
                                Ki::Number,
                                Kd::Number,
                                dt::Number,
                            intRes::Number)
    Ftipmat = Ftipmat'
    thetamatd = thetamatd'
    dthetamatd = dthetamatd'
    ddthetamatd = ddthetamatd'

    m, n = size(thetamatd)

    thetacurrent = copy(thetalist)
    dthetacurrent = copy(dthetalist)

    eint = reshape(zeros(m, 1), (m,))
    taumat = zeros(size(thetamatd))
    thetamat = zeros(size(thetamatd))

    for i = 1:n
        taulist = ComputedTorque(thetacurrent, dthetacurrent, eint, gtilde, Mtildelist, Gtildelist, Slist, thetamatd[:, i], dthetamatd[:, i], ddthetamatd[:, i], Kp, Ki, Kd)

        for j = 1:intRes
            ddthetalist = ForwardDynamics(thetacurrent, dthetacurrent, taulist, g, Ftipmat[:, i], Mlist, Glist, Slist)
            thetacurrent, dthetacurrent = EulerStep(thetacurrent, dthetacurrent, ddthetalist, dt / intRes)
        end

        taumat[:, i] = taulist
        thetamat[:, i] = thetacurrent

        eint += dt * (thetamatd[:, i] - thetacurrent)
    end

    taumat', thetamat'
end

end # module

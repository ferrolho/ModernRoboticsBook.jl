module ModernRobotics

using LinearAlgebra

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
       TestIfSE3

"""
*** BASIC HELPER FUNCTIONS ***
"""

"""
    NearZero(z)

Determines whether a scalar is small enough to be treated as zero.

# Examples
```jldoctest
julia> NearZero(-1e-7)
true
```
"""
NearZero(z::Number) = abs(z) < 1e-6

"""
    Normalize(V)

Normalizes a vector.

# Examples
```jldoctest
julia> Normalize([1 2 3])
1×3 Array{Float64,2}:
 0.267261  0.534522  0.801784
"""
Normalize(V::Array) = V / norm(V)

"""
*** CHAPTER 3: RIGID-BODY MOTIONS ***
"""

"""
    RotInv(R)

Inverts a rotation matrix.

# Examples
```jldoctest
julia> RotInv([0 0 1; 1 0 0; 0 1 0])
3×3 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
 0  1  0
 0  0  1
 1  0  0
"""
RotInv(R::Array) = R'

"""
    VecToso3(omg)

Converts a 3-vector to an so(3) representation.

# Examples
```jldoctest
julia> VecToso3([1 2 3])
3×3 Array{Int64,2}:
  0  -3   2
  3   0  -1
 -2   1   0
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
```jldoctest
julia> so3ToVec([0 -3 2; 3 0 -1; -2 1 0])
3-element Array{Int64,1}:
 1
 2
 3
"""
function so3ToVec(so3mat::Array)
    [so3mat[3, 2], so3mat[1, 3], so3mat[2, 1]]
end

"""
    AxisAng3(expc3)

Converts a 3-vector of exponential coordinates for rotation into axis-angle form.

# Examples
```jldoctest
julia> AxisAng3([1 2 3])
([0.267261 0.534522 0.801784], 3.7416573867739413)
"""
AxisAng3(expc3::Array) = Normalize(expc3), norm(expc3)

"""
    MatrixExp3(so3mat)

Computes the matrix exponential of a matrix in so(3).

# Examples
```jldoctest
julia> MatrixExp3([0 -3 2; 3 0 -1; -2 1 0])
3×3 Array{Float64,2}:
 -0.694921   0.713521  0.0892929
 -0.192007  -0.303785  0.933192 
  0.692978   0.63135   0.348107 
"""
function MatrixExp3(so3mat::Array)
    omgtheta = so3ToVec(so3mat)
    if NearZero(norm(omgtheta))
        return I
    else
        θ = AxisAng3(omgtheta)[2]
        omgmat = so3mat / θ
        return I + sin(θ) * omgmat + (1 - cos(θ)) * omgmat * omgmat
    end
end

"""
    MatrixLog3(R)

Computes the matrix logarithm of a rotation matrix.

# Examples
```jldoctest
julia> MatrixLog3([0 0 1; 1 0 0; 0 1 0])
3×3 Array{Float64,2}:
  0.0     -1.2092   1.2092
  1.2092   0.0     -1.2092
 -1.2092   1.2092   0.0   
"""
function MatrixLog3(R::Array)
    acosinput = (tr(R) - 1) / 2.0
    if acosinput >= 1
        return zeros(3, 3)
    elseif acosinput <= -1
        if not NearZero(1 + R[3, 3])
            omg = (1.0 / √(2 * (1 + R[3, 3]))) * [R[1, 3], R[2, 3], 1 + R[3, 3]]
        elseif not NearZero(1 + R[2, 2])
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
```jldoctest
julia> RpToTrans([1 0 0; 0 0 -1; 0 1 0], [1, 2, 5])
4×4 Array{Int64,2}:
 1  0   0  1
 0  0  -1  2
 0  1   0  5
 0  0   0  1
"""
RpToTrans(R::Array, p::Array) = vcat(hcat(R, p), [0 0 0 1])

"""
    TransToRp(T)

Converts a homogeneous transformation matrix into a rotation matrix and position vector.

# Examples
```jldoctest
julia> TransToRp([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
([1 0 0; 0 0 -1; 0 1 0], [0, 0, 3])
"""
TransToRp(T::Array) = T[1:3, 1:3], T[1:3, 4]

"""
    TransInv(T)

Inverts a homogeneous transformation matrix.

# Examples
```jldoctest
julia> TransInv([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
4×4 Array{Int64,2}:
 1   0  0   0
 0   0  1  -3
 0  -1  0   0
 0   0  0   1
"""
function TransInv(T::Array)
    R, p = TransToRp(T)
    vcat(hcat(R', -R' * p), [0 0 0 1])
end

"""
    VecTose3(V)

Converts a spatial velocity vector into a 4x4 matrix in se3.

# Examples
```jldoctest
julia> VecTose3([1 2 3 4 5 6])
4×4 Array{Float64,2}:
  0.0  -3.0   2.0  4.0
  3.0   0.0  -1.0  5.0
 -2.0   1.0   0.0  6.0
  0.0   0.0   0.0  0.0
"""
VecTose3(V::Array) = vcat(hcat(VecToso3(V[1:3]), V[4:6]), zeros(1, 4))

"""
    se3ToVec(se3mat)

Converts an se3 matrix into a spatial velocity vector.

# Examples
```jldoctest
julia> se3ToVec([0 -3 2 4; 3 0 -1 5; -2 1 0 6; 0 0 0 0])
1×6 Array{Int64,2}:
 1  2  3  4  5  6
"""
se3ToVec(se3mat::Array) = hcat([se3mat[3, 2] se3mat[1, 3] se3mat[2, 1]], se3mat[1:3, 4]')

"""
    Adjoint(T)

Computes the adjoint representation of a homogeneous transformation matrix.

# Examples
```jldoctest
julia> Adjoint([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
6×6 Array{Float64,2}:
 1.0  0.0   0.0  0.0  0.0   0.0
 0.0  0.0  -1.0  0.0  0.0   0.0
 0.0  1.0   0.0  0.0  0.0   0.0
 0.0  0.0   3.0  1.0  0.0   0.0
 3.0  0.0   0.0  0.0  0.0  -1.0
 0.0  0.0   0.0  0.0  1.0   0.0
"""
function Adjoint(T::Array)
    R, p = TransToRp(T)
    vcat(hcat(R, zeros(3, 3)), hcat(VecToso3(p) * R, R))
end

"""
    ScrewToAxis(q, s, h)

Takes a parametric description of a screw axis and converts it to a normalized screw axis.

# Examples
```jldoctest
julia> ScrewToAxis([3; 0; 0], [0; 0; 1], 2)
6-element Array{Int64,1}:
  0
  0
  1
  0
 -3
  2
"""
ScrewToAxis(q::Array, s::Array, h::Number) = vcat(s, q × s + h * s)

"""
    AxisAng6(expc6)

Converts a 6-vector of exponential coordinates into screw axis-angle form.

# Examples
```jldoctest
julia> AxisAng6([1, 0, 0, 1, 2, 3])
([1.0, 0.0, 0.0, 1.0, 2.0, 3.0], 1.0)
"""
function AxisAng6(expc6::Array)
    θ = norm(expc6[1:3])
    if NearZero(θ)
        θ = norm(expc6[3:6])
    end
    expc6 / θ, θ
end

"""
    MatrixExp6(se3mat)

Computes the matrix exponential of an se3 representation of exponential coordinates.

# Examples
```jldoctest
julia> MatrixExp6([0 0 0 0; 0 0 -1.57079632 2.35619449; 0 1.57079632 0 2.35619449; 0 0 0 0])
4×4 Array{Float64,2}:
 1.0  0.0         0.0        0.0       
 0.0  6.7949e-9  -1.0        1.01923e-8
 0.0  1.0         6.7949e-9  3.0       
 0.0  0.0         0.0        1.0       
"""
function MatrixExp6(se3mat::Array)
    omgtheta = so3ToVec(se3mat[1:3, 1:3])
    if NearZero(norm(omgtheta))
        return vcat(hcat(I, se3mat[1:3, 4]), [0 0 0 1])
    else
        θ = AxisAng3(omgtheta)[2]
        omgmat = se3mat[1:3, 1:3] / θ
        return vcat(hcat(MatrixExp3(se3mat[1:3, 1:3]),
                         (I * θ +
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
```jldoctest
julia> MatrixLog6([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])
4×4 Array{Float64,2}:
 0.0  0.0      0.0     0.0    
 0.0  0.0     -1.5708  2.35619
 0.0  1.5708   0.0     2.35619
 0.0  0.0      0.0     0.0    
"""
function MatrixLog6(T::Array)
    R, p = TransToRp(T)
    omgmat = MatrixLog3(R)
    if omgmat == zeros(3, 3)
        return vcat(hcat(zeros(3, 3), T[1:3, 4]), [0 0 0 0])
    else
        θ = acos((tr(R) - 1) / 2.0)
        return vcat(hcat(omgmat,
                         (I - omgmat / 2.0 +
                          (1.0 / θ - 1.0 / tan(θ / 2.0) / 2) *
                          omgmat * omgmat / θ) * T[1:3, 4]),
                    [0 0 0 0])
    end
end

"""
    ProjectToSO3(mat)

Returns a projection of mat into SO(3).

# Examples
```jldoctest
julia> ProjectToSO3([0.675 0.150  0.720; 0.370 0.771 -0.511; -0.630 0.619  0.472])
3×3 Array{Float64,2}:
  0.679011  0.148945   0.718859
  0.373207  0.773196  -0.512723
 -0.632187  0.616428   0.469421
"""
function ProjectToSO3(mat::Array)
    F  = svd(mat)
    R = F.U * F.Vt
    if det(R) < 0
        # In this case the result may be far from mat.
        R[:, F.s[3, 3]] = -R[:, F.s[3, 3]]
    end
    return R
end

"""
    ProjectToSE3(mat)

Returns a projection of mat into SE(3).

# Examples
```jldoctest
julia> ProjectToSE3([0.675 0.150 0.720 1.2; 0.370 0.771 -0.511 5.4; -0.630 0.619 0.472 3.6; 0.003 0.002 0.010 0.9])
4×4 Array{Float64,2}:
  0.679011  0.148945   0.718859  1.2
  0.373207  0.773196  -0.512723  5.4
 -0.632187  0.616428   0.469421  3.6
  0.0       0.0        0.0       1.0
"""
ProjectToSE3(mat::Array) = RpToTrans(ProjectToSO3(mat[1:3, 1:3]), mat[1:3, 4])

"""
    DistanceToSO3(mat)

Returns the Frobenius norm to describe the distance of mat from the SO(3) manifold.

# Examples
```jldoctest
julia> DistanceToSO3([1.0 0.0 0.0 ; 0.0 0.1 -0.95; 0.0 1.0 0.1])
0.08835298523536149
"""
DistanceToSO3(mat::Array) = det(mat) > 0 ? norm(mat'mat - I) : 1e+9

"""
    DistanceToSE3(mat)

Returns the Frobenius norm to describe the distance of mat from the SE(3) manifold.

# Examples
```jldoctest
julia> DistanceToSE3([1.0 0.0 0.0 1.2; 0.0 0.1 -0.95 1.5; 0.0 1.0 0.1 -0.9; 0.0 0.0 0.1 0.98])
0.13493053768513638
"""
function DistanceToSE3(mat::Array)
    matR = mat[1:3, 1:3]
    if det(matR) > 0
        norm(hcat(vcat(matR'matR, zeros(1, 3)), mat[4, :]) - I)
    else
        1e+9
    end
end

"""
    TestIfSO3(mat)

Returns true if mat is close to or on the manifold SO(3).

# Examples
```jldoctest
julia> TestIfSO3([1.0 0.0 0.0; 0.0 0.1 -0.95; 0.0 1.0 0.1])
false
"""
TestIfSO3(mat::Array) = abs(DistanceToSO3(mat)) < 1e-3

"""
    TestIfSE3(mat)

Returns true if mat is close to or on the manifold SE(3).

# Examples
```jldoctest
julia> TestIfSE3([1.0 0.0 0.0 1.2; 0.0 0.1 -0.95 1.5; 0.0 1.0 0.1 -0.9; 0.0 0.0 0.1 0.98])
false
"""
TestIfSE3(mat::Array) = abs(DistanceToSE3(mat)) < 1e-3

end # module

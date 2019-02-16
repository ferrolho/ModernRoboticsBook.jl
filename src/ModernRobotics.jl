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
       RpToTrans

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
AxisAng3(expc3::Array) = (Normalize(expc3), norm(expc3))

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
        theta = AxisAng3(omgtheta)[2]
        omgmat = so3mat / theta
        return I + sin(theta) * omgmat + (1 - cos(theta)) * omgmat * omgmat
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
        theta = acos(acosinput)
        return theta / 2.0 / sin(theta) * (R - R')
    end
end

"""
    RpToTrans(R, p)

Converts a rotation matrix and a position vector into homogeneous transformation matrix.

# Examples
```jldoctest
julia> RpToTrans([1 0 0; 0 0 -1; 0 1 0], [1 2 5])
4×4 Array{Int64,2}:
 1  0   0  1
 0  0  -1  2
 0  1   0  5
 0  0   0  1
"""
RpToTrans(R::Array, p::Array) = vcat(hcat(R, p'), [0 0 0 1])

end # module

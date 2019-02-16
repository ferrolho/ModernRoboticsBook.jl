module ModernRobotics

using LinearAlgebra

export NearZero,
       Normalize

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

end # module

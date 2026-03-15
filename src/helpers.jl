export near_zero

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

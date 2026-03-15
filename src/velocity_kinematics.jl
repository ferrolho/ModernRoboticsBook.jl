export jacobian_body, jacobian_body!, jacobian_space, jacobian_space!

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
    Jb = similar(body_screw_axes, Float64)
    jacobian_body!(Jb, body_screw_axes, joint_positions)
end

"""
    jacobian_body!(Jb, body_screw_axes, joint_positions)

In-place version of [`jacobian_body`](@ref). Writes the result into `Jb`.
"""
function jacobian_body!(
    Jb::AbstractMatrix,
    body_screw_axes::AbstractMatrix,
    joint_positions::AbstractVector,
)
    Jb .= body_screw_axes
    T = SMatrix{4,4,Float64}(LA.I)
    for i in Iterators.reverse(firstindex(joint_positions):(lastindex(joint_positions)-1))
        Bi1 = SVector{6}(@view body_screw_axes[:, i+1])
        T = T * matrix_exp6(vec_to_se3(Bi1 * -joint_positions[i+1]))
        Bi = SVector{6}(@view body_screw_axes[:, i])
        Jb[:, i] = adjoint_representation(T) * Bi
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
    Js = similar(screw_axes, Float64)
    jacobian_space!(Js, screw_axes, joint_positions)
end

"""
    jacobian_space!(Js, screw_axes, joint_positions)

In-place version of [`jacobian_space`](@ref). Writes the result into `Js`.
"""
function jacobian_space!(
    Js::AbstractMatrix,
    screw_axes::AbstractMatrix,
    joint_positions::AbstractVector,
)
    Js .= screw_axes
    T = SMatrix{4,4,Float64}(LA.I)
    for i = (firstindex(joint_positions)+1):lastindex(joint_positions)
        Si1 = SVector{6}(@view screw_axes[:, i-1])
        T = T * matrix_exp6(vec_to_se3(Si1 * joint_positions[i-1]))
        Si = SVector{6}(@view screw_axes[:, i])
        Js[:, i] = adjoint_representation(T) * Si
    end
    Js
end

export cubic_time_scaling,
    quintic_time_scaling, joint_trajectory, screw_trajectory, cartesian_trajectory

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
    traj = Vector{Matrix{Float64}}(undef, N)

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
    traj = Vector{Matrix{Float64}}(undef, N)

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

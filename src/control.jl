export computed_torque, simulate_control

"""
    computed_torque(joint_positions, joint_velocities, error_integral, gravity, link_frames, spatial_inertias, screw_axes, desired_joint_positions, desired_joint_velocities, desired_joint_accelerations, Kp, Ki, Kd)

Computes the joint control torques at a particular time instant using the computed
torque control law:

``\\tau = \\widehat{M}(\\theta)(\\ddot{\\theta}_d + K_p e + K_i \\int e + K_d \\dot{e}) + \\widehat{h}(\\theta, \\dot{\\theta})``

where ``e = \\theta_d - \\theta``, ``\\dot{e} = \\dot{\\theta}_d - \\dot{\\theta}``,
``\\widehat{M}`` is the model of the robot's mass matrix, and ``\\widehat{h}`` is the
model of centripetal, Coriolis, and gravitational forces.

!!! info "How does computed torque control work?"
    Computed torque control (also called inverse dynamics control) combines two components:

    1. **Feedforward**: uses the robot's dynamics model (``\\widehat{M}``, ``\\widehat{h}``)
       to compute the torques needed for the desired trajectory, cancelling out nonlinear
       effects (Coriolis, centripetal, gravity).
    2. **Feedback**: a PID law (``K_p e + K_i \\int e + K_d \\dot{e}``) corrects for
       tracking errors and model inaccuracies.

    Together, this linearises the closed-loop system so that each joint behaves like a
    simple second-order system, achieving precise trajectory tracking even at high speeds.

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
    mass_matrix_crba(joint_positions, link_frames, spatial_inertias, screw_axes) * (
        Kp * e +
        Ki * (error_integral + e) +
        Kd * (desired_joint_velocities - joint_velocities)
    ) + inverse_dynamics_rnea(
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
[`forward_dynamics_crba`](@ref) and numerical integration. The controller may use different
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

    error_integral = zeros(m)
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
            joint_accelerations = forward_dynamics_crba(
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

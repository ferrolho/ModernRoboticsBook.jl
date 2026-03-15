using JSON

export Robot, load_robot

const DEFAULT_GRAVITY = [0.0, 0.0, -9.81]

"""
    Robot

A data container holding the kinematic and dynamic parameters of a serial-chain robot
in the conventions of the Modern Robotics textbook.

# Fields
- `name::String`: robot name.
- `n_joints::Int`: number of actuated joints.
- `home_ee_pose::Matrix{Float64}`: end-effector pose ``M \\in SE(3)`` when all joint positions are zero (4×4).
- `screw_axes_space::Matrix{Float64}`: joint screw axes in the space (base) frame at the home configuration (6×n).
- `screw_axes_body::Matrix{Float64}`: joint screw axes in the end-effector (body) frame at the home configuration (6×n).
- `link_frames::Vector{Matrix{Float64}}`: transforms ``M_{i-1,i}`` between consecutive link frames (n+1 matrices).
- `spatial_inertias::Vector{Matrix{Float64}}`: spatial inertia matrices ``G_i`` for each link (n 6×6 matrices).
- `gravity::Vector{Float64}`: gravity vector (default `[0, 0, -9.81]`).
- `joint_names::Vector{String}`: joint names.
- `joint_types::Vector{Symbol}`: joint types (`:revolute` or `:prismatic`).
- `joint_limits::Vector{Tuple{Float64,Float64}}`: `(lower, upper)` position limits per joint.
"""
struct Robot
    name::String
    n_joints::Int
    home_ee_pose::Matrix{Float64}
    screw_axes_space::Matrix{Float64}
    screw_axes_body::Matrix{Float64}
    link_frames::Vector{Matrix{Float64}}
    spatial_inertias::Vector{Matrix{Float64}}
    gravity::Vector{Float64}
    joint_names::Vector{String}
    joint_types::Vector{Symbol}
    joint_limits::Vector{Tuple{Float64,Float64}}
end

"""
    load_robot(path::AbstractString; gravity=$DEFAULT_GRAVITY) -> Robot

Load a robot model from a JSON file.

The JSON file should contain the robot's kinematic and dynamic parameters serialized
by the companion Python converter script (`robots/convert.py`).

The `gravity` keyword argument sets the gravity vector for the robot's operating
environment. Defaults to Earth gravity `[0, 0, -9.81]`.
"""
function load_robot(path::AbstractString; gravity::AbstractVector = DEFAULT_GRAVITY)
    json = JSON.parsefile(path)

    name = json["name"]::String
    n_joints = Int(json["n_joints"])

    home_ee_pose = _json_matrix(json["home_ee_pose"])
    screw_axes_space = Matrix(_json_matrix(json["screw_axes_space"])')  # n×6 in JSON → 6×n
    screw_axes_body = Matrix(_json_matrix(json["screw_axes_body"])')
    link_frames = [_json_matrix(m) for m in json["link_frames"]]
    spatial_inertias = [_json_matrix(m) for m in json["spatial_inertias"]]

    joint_names = String[s::String for s in json["joint_names"]]
    joint_types = Symbol[Symbol(s) for s in json["joint_types"]]
    joint_limits =
        Tuple{Float64,Float64}[(Float64(l[1]), Float64(l[2])) for l in json["joint_limits"]]

    Robot(
        name,
        n_joints,
        home_ee_pose,
        screw_axes_space,
        screw_axes_body,
        link_frames,
        spatial_inertias,
        Float64.(gravity),
        joint_names,
        joint_types,
        joint_limits,
    )
end

"""
    load_robot(name::Symbol; gravity=$DEFAULT_GRAVITY) -> Robot

Load a bundled robot model by name (e.g., `load_robot(:ur5)`).

Bundled models are stored in the `robots/models/` directory of this package.
"""
function load_robot(name::Symbol; gravity::AbstractVector = DEFAULT_GRAVITY)
    path = joinpath(@__DIR__, "..", "robots", "models", "$(name).json")
    isfile(path) || error("No bundled robot model found for :$name. Expected file: $path")
    load_robot(path; gravity)
end

function _json_matrix(rows)
    reduce(vcat, [Float64.(r)' for r in rows])
end

# Convenience wrappers: kinematics

forward_kinematics_body(robot::Robot, joint_positions::AbstractVector) =
    forward_kinematics_body(
        copy(robot.home_ee_pose),
        robot.screw_axes_body,
        joint_positions,
    )

forward_kinematics_space(robot::Robot, joint_positions::AbstractVector) =
    forward_kinematics_space(
        copy(robot.home_ee_pose),
        robot.screw_axes_space,
        joint_positions,
    )

jacobian_body(robot::Robot, joint_positions::AbstractVector) =
    jacobian_body(robot.screw_axes_body, joint_positions)

jacobian_space(robot::Robot, joint_positions::AbstractVector) =
    jacobian_space(robot.screw_axes_space, joint_positions)

inverse_kinematics_body(
    robot::Robot,
    target_config::AbstractMatrix,
    initial_guess::AbstractVector,
    angular_tolerance::Number,
    linear_tolerance::Number,
) = inverse_kinematics_body(
    robot.screw_axes_body,
    robot.home_ee_pose,
    target_config,
    initial_guess,
    angular_tolerance,
    linear_tolerance,
)

inverse_kinematics_space(
    robot::Robot,
    target_config::AbstractMatrix,
    initial_guess::AbstractVector,
    angular_tolerance::Number,
    linear_tolerance::Number,
) = inverse_kinematics_space(
    robot.screw_axes_space,
    robot.home_ee_pose,
    target_config,
    initial_guess,
    angular_tolerance,
    linear_tolerance,
)

# Convenience wrappers: dynamics

function inverse_dynamics(
    robot::Robot,
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_accelerations::AbstractVector;
    gravity::AbstractVector = robot.gravity,
    tip_wrench::AbstractVector = zeros(6),
)
    inverse_dynamics(
        joint_positions,
        joint_velocities,
        joint_accelerations,
        gravity,
        tip_wrench,
        robot.link_frames,
        robot.spatial_inertias,
        robot.screw_axes_space,
    )
end

mass_matrix(robot::Robot, joint_positions::AbstractVector) = mass_matrix(
    joint_positions,
    robot.link_frames,
    robot.spatial_inertias,
    robot.screw_axes_space,
)

velocity_quadratic_forces(
    robot::Robot,
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
) = velocity_quadratic_forces(
    joint_positions,
    joint_velocities,
    robot.link_frames,
    robot.spatial_inertias,
    robot.screw_axes_space,
)

function gravity_forces(
    robot::Robot,
    joint_positions::AbstractVector;
    gravity::AbstractVector = robot.gravity,
)
    gravity_forces(
        joint_positions,
        gravity,
        robot.link_frames,
        robot.spatial_inertias,
        robot.screw_axes_space,
    )
end

function end_effector_forces(
    robot::Robot,
    joint_positions::AbstractVector,
    tip_wrench::AbstractVector,
)
    end_effector_forces(
        joint_positions,
        tip_wrench,
        robot.link_frames,
        robot.spatial_inertias,
        robot.screw_axes_space,
    )
end

function forward_dynamics(
    robot::Robot,
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_torques::AbstractVector;
    gravity::AbstractVector = robot.gravity,
    tip_wrench::AbstractVector = zeros(6),
)
    forward_dynamics(
        joint_positions,
        joint_velocities,
        joint_torques,
        gravity,
        tip_wrench,
        robot.link_frames,
        robot.spatial_inertias,
        robot.screw_axes_space,
    )
end

# Convenience wrappers: trajectory dynamics

function inverse_dynamics_trajectory(
    robot::Robot,
    joint_position_traj::AbstractMatrix,
    joint_velocity_traj::AbstractMatrix,
    joint_acceleration_traj::AbstractMatrix;
    gravity::AbstractVector = robot.gravity,
    tip_wrench_traj::AbstractMatrix = zeros(size(joint_position_traj, 1), 6),
)
    inverse_dynamics_trajectory(
        joint_position_traj,
        joint_velocity_traj,
        joint_acceleration_traj,
        gravity,
        tip_wrench_traj,
        robot.link_frames,
        robot.spatial_inertias,
        robot.screw_axes_space,
    )
end

function forward_dynamics_trajectory(
    robot::Robot,
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_torque_traj::AbstractMatrix,
    timestep::Number,
    integration_resolution::Number;
    gravity::AbstractVector = robot.gravity,
    tip_wrench_traj::AbstractMatrix = zeros(size(joint_torque_traj, 1), 6),
)
    forward_dynamics_trajectory(
        joint_positions,
        joint_velocities,
        joint_torque_traj,
        gravity,
        tip_wrench_traj,
        robot.link_frames,
        robot.spatial_inertias,
        robot.screw_axes_space,
        timestep,
        integration_resolution,
    )
end

# Convenience wrappers: control

function computed_torque(
    robot::Robot,
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    error_integral::AbstractVector,
    desired_joint_positions::AbstractVector,
    desired_joint_velocities::AbstractVector,
    desired_joint_accelerations::AbstractVector,
    Kp::Number,
    Ki::Number,
    Kd::Number;
    gravity::AbstractVector = robot.gravity,
)
    computed_torque(
        joint_positions,
        joint_velocities,
        error_integral,
        gravity,
        robot.link_frames,
        robot.spatial_inertias,
        robot.screw_axes_space,
        desired_joint_positions,
        desired_joint_velocities,
        desired_joint_accelerations,
        Kp,
        Ki,
        Kd,
    )
end

function simulate_control(
    robot::Robot,
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    tip_wrench_traj::AbstractMatrix,
    desired_joint_position_traj::AbstractMatrix,
    desired_joint_velocity_traj::AbstractMatrix,
    desired_joint_acceleration_traj::AbstractMatrix,
    estimated_robot::Robot,
    Kp::Number,
    Ki::Number,
    Kd::Number,
    timestep::Number,
    integration_resolution::Number;
    gravity::AbstractVector = robot.gravity,
    estimated_gravity::AbstractVector = estimated_robot.gravity,
)
    simulate_control(
        joint_positions,
        joint_velocities,
        gravity,
        tip_wrench_traj,
        robot.link_frames,
        robot.spatial_inertias,
        robot.screw_axes_space,
        desired_joint_position_traj,
        desired_joint_velocity_traj,
        desired_joint_acceleration_traj,
        estimated_gravity,
        estimated_robot.link_frames,
        estimated_robot.spatial_inertias,
        Kp,
        Ki,
        Kd,
        timestep,
        integration_resolution,
    )
end

#!/usr/bin/env python3
"""Convert URDF robot descriptions to Modern Robotics parameter JSON files.

This script extracts the kinematic and dynamic parameters needed by
ModernRoboticsBook.jl from standard URDF files. It produces JSON files
that can be loaded with `load_robot()` in Julia.

Dependencies:
    pip install yourdfpy numpy robot-descriptions

Usage:
    # Convert a bundled robot by name
    python convert.py ur5

    # Convert a specific URDF file
    python convert.py --urdf path/to/robot.urdf --name my_robot

    # List available bundled robots
    python convert.py --list
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

# Available robots from robot_descriptions.py
ROBOT_REGISTRY = {
    "ur3": "ur3_description",
    "ur3e": "ur3e_description",
    "ur5": "ur5_description",
    "ur5e": "ur5e_description",
    "ur10": "ur10_description",
    "ur10e": "ur10e_description",
    "panda": "panda_description",
}


def adjoint(T):
    """Compute the 6x6 adjoint representation of a 4x4 SE(3) transform."""
    R = T[:3, :3]
    p = T[:3, 3]
    p_hat = np.array([
        [0, -p[2], p[1]],
        [p[2], 0, -p[0]],
        [-p[1], p[0], 0],
    ])
    Ad = np.zeros((6, 6))
    Ad[:3, :3] = R
    Ad[3:, :3] = p_hat @ R
    Ad[3:, 3:] = R
    return Ad


def trans_inv(T):
    """Compute the inverse of an SE(3) transform."""
    R = T[:3, :3]
    p = T[:3, 3]
    T_inv = np.eye(4)
    T_inv[:3, :3] = R.T
    T_inv[:3, 3] = -R.T @ p
    return T_inv


def urdf_to_mr(urdf_path):
    """Extract Modern Robotics parameters from a URDF file.

    Uses yourdfpy for URDF parsing. Joint/inertial origins are already
    parsed as 4x4 SE(3) transforms by yourdfpy.

    Returns a dict with keys: name, n_joints, home_config, screw_axes_space,
    screw_axes_body, link_frames, spatial_inertias, joint_names, joint_types,
    joint_limits.
    """
    import yourdfpy

    robot = yourdfpy.URDF.load(str(urdf_path))
    urdf = robot.robot

    # Build parent-child maps
    child_to_joint = {j.child: j for j in urdf.joints}
    parent_to_joints = {}
    for j in urdf.joints:
        parent_to_joints.setdefault(j.parent, []).append(j)

    # Find the root link (not a child of any joint)
    all_children = {j.child for j in urdf.joints}
    all_parents = {j.parent for j in urdf.joints}
    root_links = all_parents - all_children
    if len(root_links) != 1:
        raise ValueError(f"Expected 1 root link, found {len(root_links)}: {root_links}")
    root_link = root_links.pop()

    # Walk the kinematic tree from root, collecting actuated joints in order
    actuated_joints = []

    def walk(link_name):
        for j in parent_to_joints.get(link_name, []):
            if j.type in ("revolute", "prismatic", "continuous"):
                actuated_joints.append(j)
            walk(j.child)

    walk(root_link)

    n = len(actuated_joints)
    if n == 0:
        raise ValueError("No actuated joints found in URDF")

    # Compute FK at zero configuration: cumulative transform from base to each link
    # yourdfpy already gives us origin as a 4x4 SE(3) matrix
    link_transforms = {root_link: np.eye(4)}

    def compute_fk(link_name):
        for j in parent_to_joints.get(link_name, []):
            T_parent = link_transforms[link_name]
            T_joint = j.origin  # already a 4x4 numpy array from yourdfpy
            link_transforms[j.child] = T_parent @ T_joint
            compute_fk(j.child)

    compute_fk(root_link)

    # Find end-effector link: walk past fixed joints after the last actuated joint
    ee_link = actuated_joints[-1].child
    current = ee_link
    while current in parent_to_joints:
        fixed_children = [j for j in parent_to_joints[current] if j.type == "fixed"]
        if len(fixed_children) == 1:
            current = fixed_children[0].child
        else:
            break
    ee_link = current

    # M: end-effector home configuration
    M = link_transforms[ee_link]

    # Slist: space-frame screw axes
    Slist = np.zeros((6, n))
    for i, joint in enumerate(actuated_joints):
        # Joint frame in the space frame
        T_parent = link_transforms[joint.parent]
        T_joint_frame = T_parent @ joint.origin
        R = T_joint_frame[:3, :3]
        p = T_joint_frame[:3, 3]

        # Joint axis in space frame
        axis = np.array(joint.axis if joint.axis is not None else [0, 0, 1], dtype=float)
        w = R @ axis
        w = w / np.linalg.norm(w)

        if joint.type in ("revolute", "continuous"):
            v = -np.cross(w, p)
            Slist[:, i] = np.concatenate([w, v])
        elif joint.type == "prismatic":
            Slist[:, i] = np.concatenate([np.zeros(3), w])

    # Blist: body-frame screw axes
    Blist = adjoint(trans_inv(M)) @ Slist

    # Link inertial properties
    link_map = {link.name: link for link in urdf.links}

    # Mlist: relative transforms between consecutive link CoM frames
    com_transforms = []
    for joint in actuated_joints:
        child_link = link_map[joint.child]
        T_0_link = link_transforms[joint.child]

        if child_link.inertial is not None and child_link.inertial.origin is not None:
            T_inertial = child_link.inertial.origin  # already 4x4 from yourdfpy
        else:
            T_inertial = np.eye(4)

        T_0_com = T_0_link @ T_inertial
        com_transforms.append(T_0_com)

    Mlist = []
    Mlist.append(com_transforms[0])
    for i in range(1, n):
        Mlist.append(trans_inv(com_transforms[i - 1]) @ com_transforms[i])
    Mlist.append(trans_inv(com_transforms[-1]) @ M)

    # Glist: spatial inertia matrices
    Glist = []
    for joint in actuated_joints:
        child_link = link_map[joint.child]
        G = np.zeros((6, 6))

        if child_link.inertial is not None:
            inertial = child_link.inertial
            mass = inertial.mass if inertial.mass is not None else 0.0

            if inertial.inertia is not None:
                # yourdfpy gives inertia as a 3x3 numpy array
                G[:3, :3] = inertial.inertia

                # Rotate inertia if inertial origin has non-identity rotation
                if inertial.origin is not None:
                    R_inertial = inertial.origin[:3, :3]
                    if not np.allclose(R_inertial, np.eye(3)):
                        G[:3, :3] = R_inertial @ G[:3, :3] @ R_inertial.T

            G[3, 3] = mass
            G[4, 4] = mass
            G[5, 5] = mass

        Glist.append(G)

    # Joint metadata
    joint_names = [j.name for j in actuated_joints]
    joint_types = []
    joint_limits = []
    for j in actuated_joints:
        jtype = "revolute" if j.type in ("revolute", "continuous") else "prismatic"
        joint_types.append(jtype)

        if j.limit is not None:
            lower = j.limit.lower if j.limit.lower is not None else -6.283185307
            upper = j.limit.upper if j.limit.upper is not None else 6.283185307
        else:
            lower, upper = -6.283185307, 6.283185307
        joint_limits.append([lower, upper])

    return {
        "name": urdf.name or Path(urdf_path).stem,
        "n_joints": n,
        "home_config": M.tolist(),
        "screw_axes_space": Slist.T.tolist(),  # n x 6 (each row is a screw axis)
        "screw_axes_body": Blist.T.tolist(),
        "link_frames": [m.tolist() for m in Mlist],
        "spatial_inertias": [g.tolist() for g in Glist],
        "joint_names": joint_names,
        "joint_types": joint_types,
        "joint_limits": joint_limits,
    }


def get_urdf_path(robot_name):
    """Get URDF path for a named robot from robot_descriptions.py."""
    if robot_name not in ROBOT_REGISTRY:
        available = ", ".join(sorted(ROBOT_REGISTRY.keys()))
        raise ValueError(
            f"Unknown robot '{robot_name}'. Available: {available}\n"
            f"Use --urdf to specify a custom URDF file."
        )

    module_name = ROBOT_REGISTRY[robot_name]
    import importlib
    mod = importlib.import_module(f"robot_descriptions.{module_name}")
    return mod.URDF_PATH


def main():
    parser = argparse.ArgumentParser(
        description="Convert URDF to Modern Robotics parameter JSON files."
    )
    parser.add_argument(
        "robot", nargs="?", help="Robot name from the registry (e.g., ur5, panda)"
    )
    parser.add_argument(
        "--urdf", help="Path to a custom URDF file"
    )
    parser.add_argument(
        "--name", help="Robot name for the output (defaults to URDF filename)"
    )
    parser.add_argument(
        "--output", "-o", help="Output JSON file path (defaults to robots/models/<name>.json)"
    )
    parser.add_argument(
        "--list", action="store_true", help="List available robots"
    )

    args = parser.parse_args()

    if args.list:
        print("Available robots:")
        for name in sorted(ROBOT_REGISTRY.keys()):
            print(f"  {name}")
        return

    if args.urdf:
        urdf_path = args.urdf
        robot_name = args.name or Path(urdf_path).stem
    elif args.robot:
        urdf_path = get_urdf_path(args.robot)
        robot_name = args.robot
    else:
        parser.error("Specify a robot name or --urdf path. Use --list to see available robots.")
        return

    print(f"Converting {robot_name} from {urdf_path}...")
    result = urdf_to_mr(urdf_path)

    if args.name:
        result["name"] = args.name

    output_dir = Path(__file__).parent / "models"
    output_dir.mkdir(exist_ok=True)
    output_path = args.output or str(output_dir / f"{robot_name}.json")

    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"Written to {output_path}")
    print(f"  Joints: {result['n_joints']}")
    print(f"  Joint names: {result['joint_names']}")


if __name__ == "__main__":
    main()

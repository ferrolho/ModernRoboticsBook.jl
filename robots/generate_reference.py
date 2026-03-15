#!/usr/bin/env python3
"""Generate reference values from Pinocchio for cross-validating ModernRoboticsBook.jl.

This script loads robot URDFs via Pinocchio and computes FK, Jacobians, inverse
dynamics, mass matrix, and gravity forces at several test configurations. Results
are saved as JSON for comparison in Julia tests.

Dependencies:
    pip install pin robot-descriptions numpy

Usage:
    python generate_reference.py
"""

import json
from pathlib import Path

import numpy as np
import pinocchio as pin
from robot_descriptions import ur5_description


def compute_reference(model, data, ee_frame_id, q, v, a, gravity):
    """Compute all reference values at a given configuration."""
    # Forward kinematics
    pin.framesForwardKinematics(model, data, q)
    ee_pose = data.oMf[ee_frame_id]
    T = np.eye(4)
    T[:3, :3] = ee_pose.rotation
    T[:3, 3] = ee_pose.translation

    # Jacobians
    pin.computeJointJacobians(model, data, q)
    # LOCAL = body frame Jacobian, WORLD = space frame Jacobian
    J_body = pin.getFrameJacobian(model, data, ee_frame_id, pin.LOCAL)
    J_space = pin.getFrameJacobian(model, data, ee_frame_id, pin.WORLD)

    # Mass matrix
    M = pin.crba(model, data, q)
    # crba only fills the upper triangle, symmetrise
    M = np.triu(M) + np.triu(M, 1).T

    # Gravity forces
    g = pin.computeGeneralizedGravity(model, data, q)

    # Inverse dynamics: tau = M*a + c(q,v) + g(q)
    tau = pin.rnea(model, data, q, v, a)

    # Coriolis/centrifugal (bias - gravity)
    bias = pin.rnea(model, data, q, v, np.zeros(model.nv))
    coriolis = bias - g

    return {
        "ee_pose": T.tolist(),
        "jacobian_body": J_body.tolist(),
        "jacobian_space": J_space.tolist(),
        "mass_matrix": M.tolist(),
        "gravity_forces": g.tolist(),
        "inverse_dynamics": tau.tolist(),
        "coriolis_forces": coriolis.tolist(),
    }


def main():
    # Load UR5
    model = pin.buildModelFromUrdf(ur5_description.URDF_PATH)
    data = model.createData()

    # Use ee_link to match what our URDF converter produces (standard ROS convention)
    ee_frame_id = model.getFrameId("ee_link")
    print(f"End-effector frame: ee_link (id={ee_frame_id})")
    print(f"Joints: {model.nq}")
    print(f"Gravity: {model.gravity.linear}")

    # Test configurations
    configs = {
        "zero": {
            "q": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "v": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "a": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        },
        "shoulder_lift": {
            "q": [0.0, -1.5707963267948966, 0.0, 0.0, 0.0, 0.0],
            "v": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "a": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        },
        "general": {
            "q": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            "v": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            "a": [1.0, 1.5, 2.0, 0.5, 0.3, 0.1],
        },
    }

    result = {
        "robot": "ur5",
        "source": "pinocchio",
        "ee_frame": "ee_link",
        "gravity": model.gravity.linear.tolist(),
        "n_joints": model.nq,
        "configurations": {},
    }

    for name, cfg in configs.items():
        q = np.array(cfg["q"])
        v = np.array(cfg["v"])
        a = np.array(cfg["a"])
        ref = compute_reference(model, data, ee_frame_id, q, v, a, model.gravity)
        result["configurations"][name] = {
            "q": cfg["q"],
            "v": cfg["v"],
            "a": cfg["a"],
            **ref,
        }
        print(f"\n--- {name} ---")
        print(f"  q = {cfg['q']}")
        print(f"  EE position: {ref['ee_pose'][0][3]:.6f}, {ref['ee_pose'][1][3]:.6f}, {ref['ee_pose'][2][3]:.6f}")

    output_path = Path(__file__).parent / "models" / "ur5_pinocchio_reference.json"
    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"\nWritten to {output_path}")


if __name__ == "__main__":
    main()

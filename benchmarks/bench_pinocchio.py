#!/usr/bin/env python3
"""Benchmark Pinocchio on UR5 for comparison with ModernRoboticsBook.jl."""

import timeit
import numpy as np
import pinocchio as pin
from robot_descriptions import ur5_description

model = pin.buildModelFromUrdf(ur5_description.URDF_PATH)
data = model.createData()

ee_frame_id = model.getFrameId("ee_link")

q = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
v = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
a = np.array([1.0, 1.5, 2.0, 0.5, 0.3, 0.1])

N = 100_000


def bench(label, func):
    # Warmup
    func()
    total = timeit.timeit(func, number=N)
    us = total / N * 1e6
    print(f"{label:<25s} {us:8.2f} μs")


print(f"=== Pinocchio {pin.__version__} — UR5 (6-DOF) ===\n")


def fk_space():
    pin.framesForwardKinematics(model, data, q)
    return data.oMf[ee_frame_id]


bench("FK (frames):         ", fk_space)


def jacobian_space():
    pin.computeJointJacobians(model, data, q)
    pin.framesForwardKinematics(model, data, q)
    return pin.getFrameJacobian(model, data, ee_frame_id, pin.WORLD)


bench("Jacobian (space):    ", jacobian_space)


def jacobian_body():
    pin.computeJointJacobians(model, data, q)
    pin.framesForwardKinematics(model, data, q)
    return pin.getFrameJacobian(model, data, ee_frame_id, pin.LOCAL)


bench("Jacobian (body):     ", jacobian_body)


def inverse_dynamics():
    return pin.rnea(model, data, q, v, a)


bench("Inverse dynamics:    ", inverse_dynamics)


def mass_matrix():
    M = pin.crba(model, data, q)
    return M


bench("Mass matrix:         ", mass_matrix)


def gravity_forces():
    return pin.computeGeneralizedGravity(model, data, q)


bench("Gravity forces:      ", gravity_forces)


def forward_dynamics():
    return pin.aba(model, data, q, v, np.zeros(6))


bench("Forward dynamics:    ", forward_dynamics)

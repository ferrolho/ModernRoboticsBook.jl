var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#ModernRobotics.Adjoint-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.Adjoint",
    "category": "method",
    "text": "Adjoint(T)\n\nComputes the adjoint representation of a homogeneous transformation matrix.\n\nExamples\n\njulia> Adjoint([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])\n6×6 Array{Float64,2}:\n 1.0  0.0   0.0  0.0  0.0   0.0\n 0.0  0.0  -1.0  0.0  0.0   0.0\n 0.0  1.0   0.0  0.0  0.0   0.0\n 0.0  0.0   3.0  1.0  0.0   0.0\n 3.0  0.0   0.0  0.0  0.0  -1.0\n 0.0  0.0   0.0  0.0  1.0   0.0\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.AxisAng3-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.AxisAng3",
    "category": "method",
    "text": "AxisAng3(expc3)\n\nConverts a 3-vector of exponential coordinates for rotation into axis-angle form.\n\nExamples\n\njulia> AxisAng3([1, 2, 3])\n([0.267261 0.534522 0.801784], 3.7416573867739413)\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.AxisAng6-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.AxisAng6",
    "category": "method",
    "text": "AxisAng6(expc6)\n\nConverts a 6-vector of exponential coordinates into screw axis-angle form.\n\nExamples\n\njulia> AxisAng6([1, 0, 0, 1, 2, 3])\n([1.0, 0.0, 0.0, 1.0, 2.0, 3.0], 1.0)\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.ComputedTorque-Tuple{Array,Array,Array,Array,Array,Array,AbstractArray{T,2} where T,Array,Array,Array,Number,Number,Number}",
    "page": "Home",
    "title": "ModernRobotics.ComputedTorque",
    "category": "method",
    "text": "ComputedTorque(thetalist, dthetalist, eint, g, Mlist, Glist, Slist, thetalistd, dthetalistd, ddthetalistd, Kp, Ki, Kd)\n\nComputes the joint control torques at a particular time instant.\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.CubicTimeScaling-Tuple{Number,Number}",
    "page": "Home",
    "title": "ModernRobotics.CubicTimeScaling",
    "category": "method",
    "text": "CubicTimeScaling(Tf, t)\n\nComputes s(t) for a cubic time scaling.\n\nExamples\n\njulia> CubicTimeScaling(2, 0.6)\n0.21600000000000003\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.DistanceToSE3-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.DistanceToSE3",
    "category": "method",
    "text": "DistanceToSE3(mat)\n\nReturns the Frobenius norm to describe the distance of mat from the SE(3) manifold.\n\nExamples\n\njulia> DistanceToSE3([1.0 0.0 0.0 1.2; 0.0 0.1 -0.95 1.5; 0.0 1.0 0.1 -0.9; 0.0 0.0 0.1 0.98])\n0.13493053768513638\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.DistanceToSO3-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.DistanceToSO3",
    "category": "method",
    "text": "DistanceToSO3(mat)\n\nReturns the Frobenius norm to describe the distance of mat from the SO(3) manifold.\n\nExamples\n\njulia> DistanceToSO3([1.0 0.0 0.0 ; 0.0 0.1 -0.95; 0.0 1.0 0.1])\n0.08835298523536149\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.EndEffectorForces-Tuple{Array,Array,Array,Array,AbstractArray{T,2} where T}",
    "page": "Home",
    "title": "ModernRobotics.EndEffectorForces",
    "category": "method",
    "text": "EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist)\n\nComputes the joint forces/torques an open chain robot requires only to create the end-effector force Ftip.\n\nExamples\n\njulia> EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist)\n3-element Array{Float64,1}:\n 1.4095460782639782\n 1.8577149723180628\n 1.392409          \n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.EulerStep-Tuple{Array,Array,Array,Number}",
    "page": "Home",
    "title": "ModernRobotics.EulerStep",
    "category": "method",
    "text": "EulerStep(thetalist, dthetalist, ddthetalist, dt)\n\nCompute the joint angles and velocities at the next timestep using first order Euler integration.\n\nExamples\n\njulia> EulerStep(thetalist, dthetalist, ddthetalist, dt)\n([0.11, 0.12, 0.13], [0.3, 0.35, 0.4])\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.FKinBody-Tuple{AbstractArray{T,2} where T,AbstractArray{T,2} where T,Array}",
    "page": "Home",
    "title": "ModernRobotics.FKinBody",
    "category": "method",
    "text": "FKinBody(M, Blist, thetalist)\n\nComputes forward kinematics in the body frame for an open chain robot.\n\nExamples\n\njulia> FKinBody(M, Blist, thetalist)\n4×4 Array{Float64,2}:\n -1.14424e-17  1.0           0.0  -5.0    \n  1.0          1.14424e-17   0.0   4.0    \n  0.0          0.0          -1.0   1.68584\n  0.0          0.0           0.0   1.0    \n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.FKinSpace-Tuple{AbstractArray{T,2} where T,AbstractArray{T,2} where T,Array}",
    "page": "Home",
    "title": "ModernRobotics.FKinSpace",
    "category": "method",
    "text": "FKinSpace(M, Slist, thetalist)\n\nComputes forward kinematics in the space frame for an open chain robot.\n\nExamples\n\njulia> FKinSpace(M, Slist, thetalist)\n4×4 Array{Float64,2}:\n -1.14424e-17  1.0           0.0  -5.0    \n  1.0          1.14424e-17   0.0   4.0    \n  0.0          0.0          -1.0   1.68584\n  0.0          0.0           0.0   1.0    \n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.ForwardDynamics-Tuple{Array,Array,Array,Array,Array,Array,Array,AbstractArray{T,2} where T}",
    "page": "Home",
    "title": "ModernRobotics.ForwardDynamics",
    "category": "method",
    "text": "ForwardDynamics(thetalist, dthetalist, taulist, g, Ftip, Mlist, Glist, Slist)\n\nComputes forward dynamics in the space frame for an open chain robot.\n\nExamples\n\njulia> ForwardDynamics(thetalist, dthetalist, taulist, g, Ftip, Mlist, Glist, Slist)\n3-element Array{Float64,1}:\n  -0.9739290670855626\n  25.584667840340558 \n -32.91499212478149  \n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.ForwardDynamicsTrajectory-Tuple{Array,Array,AbstractArray{T,2} where T,Array,Array,Array,Array,AbstractArray{T,2} where T,Number,Number}",
    "page": "Home",
    "title": "ModernRobotics.ForwardDynamicsTrajectory",
    "category": "method",
    "text": "ForwardDynamicsTrajectory(thetalist, dthetalist, taumat, g, Ftipmat, Mlist, Glist, Slist, dt, intRes)\n\nSimulates the motion of a serial chain given an open-loop history of joint forces/torques.\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.GravityForces-Tuple{Array,Array,Array,Array,AbstractArray{T,2} where T}",
    "page": "Home",
    "title": "ModernRobotics.GravityForces",
    "category": "method",
    "text": "GravityForces(thetalist, g, Mlist, Glist, Slist)\n\nComputes the joint forces/torques an open chain robot requires to overcome gravity at its configuration.\n\nExamples\n\njulia> GravityForces(thetalist, g, Mlist, Glist, Slist)\n3-element Array{Float64,1}:\n  28.40331261821983  \n -37.64094817177068  \n  -5.4415891999683605\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.IKinBody-Tuple{AbstractArray{T,2} where T,AbstractArray{T,2} where T,AbstractArray{T,2} where T,Array,Number,Number}",
    "page": "Home",
    "title": "ModernRobotics.IKinBody",
    "category": "method",
    "text": "IKinBody(Blist, M, T, thetalist0, eomg, ev)\n\nComputes inverse kinematics in the body frame for an open chain robot.\n\nExamples\n\njulia> IKinBody(Blist, M, T, thetalist0, eomg, ev)\n([1.57074; 2.99967; 3.14154], true)\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.IKinSpace-Tuple{AbstractArray{T,2} where T,AbstractArray{T,2} where T,AbstractArray{T,2} where T,Array,Number,Number}",
    "page": "Home",
    "title": "ModernRobotics.IKinSpace",
    "category": "method",
    "text": "IKinSpace(Slist, M, T, thetalist0, eomg, ev)\n\nComputes inverse kinematics in the space frame for an open chain robot.\n\nExamples\n\njulia> IKinSpace(Slist, M, T, thetalist0, eomg, ev)\n([1.57074; 2.99966; 3.14153], true)\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.InverseDynamics-Tuple{Array,Array,Array,Array,Array,Array,Array,AbstractArray{T,2} where T}",
    "page": "Home",
    "title": "ModernRobotics.InverseDynamics",
    "category": "method",
    "text": "InverseDynamics(thetalist, dthetalist, ddthetalist, g, Ftip, Mlist, Glist, Slist)\n\nComputes inverse dynamics in the space frame for an open chain robot.\n\nExamples\n\njulia> InverseDynamics(thetalist, dthetalist, ddthetalist, g, Ftip, Mlist, Glist, Slist)\n3-element Array{Float64,1}:\n  74.69616155287451 \n -33.06766015851458 \n  -3.230573137901424\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.InverseDynamicsTrajectory-Tuple{Array,Array,Array,Array,Array,Array,Array,AbstractArray{T,2} where T}",
    "page": "Home",
    "title": "ModernRobotics.InverseDynamicsTrajectory",
    "category": "method",
    "text": "InverseDynamicsTrajectory(thetamat, dthetamat, ddthetamat, g, Ftipmat, Mlist, Glist, Slist)\n\nCalculates the joint forces/torques required to move the serial chain along the given trajectory using inverse dynamics.\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.JacobianBody-Tuple{AbstractArray{T,2} where T,Array}",
    "page": "Home",
    "title": "ModernRobotics.JacobianBody",
    "category": "method",
    "text": "JacobianBody(Blist, thetalist)\n\nComputes the body Jacobian for an open chain robot.\n\nExamples\n\njulia> JacobianBody(Blist, thetalist)\n6×4 Array{Float64,2}:\n -0.0452841  0.995004    0.0       1.0\n  0.743593   0.0930486   0.362358  0.0\n -0.667097   0.0361754  -0.932039  0.0\n  2.32586    1.66809     0.564108  0.2\n -1.44321    2.94561     1.43307   0.3\n -2.0664     1.82882    -1.58869   0.4\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.JacobianSpace-Tuple{AbstractArray{T,2} where T,Array}",
    "page": "Home",
    "title": "ModernRobotics.JacobianSpace",
    "category": "method",
    "text": "JacobianSpace(Slist, thetalist)\n\nComputes the space Jacobian for an open chain robot.\n\nExamples\n\njulia> JacobianSpace(Slist, thetalist)\n6×4 Array{Float64,2}:\n 0.0  0.980067  -0.0901156   0.957494 \n 0.0  0.198669   0.444554    0.284876 \n 1.0  0.0        0.891207   -0.0452841\n 0.0  1.95219   -2.21635    -0.511615 \n 0.2  0.436541  -2.43713     2.77536  \n 0.2  2.96027    3.23573     2.22512  \n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.JointTrajectory-Tuple{Array,Array,Number,Integer,Integer}",
    "page": "Home",
    "title": "ModernRobotics.JointTrajectory",
    "category": "method",
    "text": "JointTrajectory(thetastart, thetaend, Tf, N, method)\n\nComputes a straight-line trajectory in joint space.\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.MassMatrix-Tuple{Array,Array,Array,AbstractArray{T,2} where T}",
    "page": "Home",
    "title": "ModernRobotics.MassMatrix",
    "category": "method",
    "text": "MassMatrix(thetalist, Mlist, Glist, Slist)\n\nComputes the mass matrix of an open chain robot based on the given configuration.\n\nExamples\n\njulia> MassMatrix(thetalist, Mlist, Glist, Slist)\n3×3 Array{Float64,2}:\n 22.5433      -0.307147  -0.00718426\n -0.307147     1.96851    0.432157  \n -0.00718426   0.432157   0.191631  \n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.MatrixExp3-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.MatrixExp3",
    "category": "method",
    "text": "MatrixExp3(so3mat)\n\nComputes the matrix exponential of a matrix in so(3).\n\nExamples\n\njulia> MatrixExp3([0 -3 2; 3 0 -1; -2 1 0])\n3×3 Array{Float64,2}:\n -0.694921   0.713521  0.0892929\n -0.192007  -0.303785  0.933192 \n  0.692978   0.63135   0.348107 \n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.MatrixExp6-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.MatrixExp6",
    "category": "method",
    "text": "MatrixExp6(se3mat)\n\nComputes the matrix exponential of an se3 representation of exponential coordinates.\n\nExamples\n\njulia> MatrixExp6([0 0 0 0; 0 0 -1.57079632 2.35619449; 0 1.57079632 0 2.35619449; 0 0 0 0])\n4×4 Array{Float64,2}:\n 1.0  0.0         0.0        0.0       \n 0.0  6.7949e-9  -1.0        1.01923e-8\n 0.0  1.0         6.7949e-9  3.0       \n 0.0  0.0         0.0        1.0       \n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.MatrixLog3-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.MatrixLog3",
    "category": "method",
    "text": "MatrixLog3(R)\n\nComputes the matrix logarithm of a rotation matrix.\n\nExamples\n\njulia> MatrixLog3([0 0 1; 1 0 0; 0 1 0])\n3×3 Array{Float64,2}:\n  0.0     -1.2092   1.2092\n  1.2092   0.0     -1.2092\n -1.2092   1.2092   0.0   \n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.MatrixLog6-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.MatrixLog6",
    "category": "method",
    "text": "MatrixLog6(T)\n\nComputes the matrix logarithm of a homogeneous transformation matrix.\n\nExamples\n\njulia> MatrixLog6([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])\n4×4 Array{Float64,2}:\n 0.0  0.0      0.0     0.0    \n 0.0  0.0     -1.5708  2.35619\n 0.0  1.5708   0.0     2.35619\n 0.0  0.0      0.0     0.0    \n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.NearZero-Tuple{Number}",
    "page": "Home",
    "title": "ModernRobotics.NearZero",
    "category": "method",
    "text": "NearZero(z)\n\nDetermines whether a scalar is small enough to be treated as zero.\n\nExamples\n\njulia> NearZero(-1e-7)\ntrue\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.Normalize-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.Normalize",
    "category": "method",
    "text": "Normalize(V)\n\nNormalizes a vector.\n\nExamples\n\njulia> Normalize([1, 2, 3])\n3-element Array{Float64,1}:\n 0.2672612419124244\n 0.5345224838248488\n 0.8017837257372732\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.ProjectToSE3-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.ProjectToSE3",
    "category": "method",
    "text": "ProjectToSE3(mat)\n\nReturns a projection of mat into SE(3).\n\nExamples\n\njulia> ProjectToSE3([0.675 0.150 0.720 1.2; 0.370 0.771 -0.511 5.4; -0.630 0.619 0.472 3.6; 0.003 0.002 0.010 0.9])\n4×4 Array{Float64,2}:\n  0.679011  0.148945   0.718859  1.2\n  0.373207  0.773196  -0.512723  5.4\n -0.632187  0.616428   0.469421  3.6\n  0.0       0.0        0.0       1.0\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.ProjectToSO3-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.ProjectToSO3",
    "category": "method",
    "text": "ProjectToSO3(mat)\n\nReturns a projection of mat into SO(3).\n\nExamples\n\njulia> ProjectToSO3([0.675 0.150  0.720; 0.370 0.771 -0.511; -0.630 0.619  0.472])\n3×3 Array{Float64,2}:\n  0.679011  0.148945   0.718859\n  0.373207  0.773196  -0.512723\n -0.632187  0.616428   0.469421\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.QuinticTimeScaling-Tuple{Number,Number}",
    "page": "Home",
    "title": "ModernRobotics.QuinticTimeScaling",
    "category": "method",
    "text": "QuinticTimeScaling(Tf, t)\n\nComputes s(t) for a quintic time scaling.\n\nExamples\n\njulia> QuinticTimeScaling(2, 0.6)\n0.16308\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.RotInv-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.RotInv",
    "category": "method",
    "text": "RotInv(R)\n\nInverts a rotation matrix.\n\nExamples\n\njulia> RotInv([0 0 1; 1 0 0; 0 1 0])\n3×3 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:\n 0  1  0\n 0  0  1\n 1  0  0\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.RpToTrans-Tuple{Array,Array}",
    "page": "Home",
    "title": "ModernRobotics.RpToTrans",
    "category": "method",
    "text": "RpToTrans(R, p)\n\nConverts a rotation matrix and a position vector into homogeneous transformation matrix.\n\nExamples\n\njulia> RpToTrans([1 0 0; 0 0 -1; 0 1 0], [1, 2, 5])\n4×4 Array{Int64,2}:\n 1  0   0  1\n 0  0  -1  2\n 0  1   0  5\n 0  0   0  1\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.ScrewToAxis-Tuple{Array,Array,Number}",
    "page": "Home",
    "title": "ModernRobotics.ScrewToAxis",
    "category": "method",
    "text": "ScrewToAxis(q, s, h)\n\nTakes a parametric description of a screw axis and converts it to a normalized screw axis.\n\nExamples\n\njulia> ScrewToAxis([3; 0; 0], [0; 0; 1], 2)\n6-element Array{Int64,1}:\n  0\n  0\n  1\n  0\n -3\n  2\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.SimulateControl-Tuple{Array,Array,Array,Array,Array,Array,AbstractArray{T,2} where T,Array,Array,Array,Array,Array,Array,Number,Number,Number,Number,Number}",
    "page": "Home",
    "title": "ModernRobotics.SimulateControl",
    "category": "method",
    "text": "SimulateControl(thetalist, dthetalist, g, Ftipmat, Mlist, Glist,\n                Slist, thetamatd, dthetamatd, ddthetamatd, gtilde,\n                Mtildelist, Gtildelist, Kp, Ki, Kd, dt, intRes)\n\nSimulates the computed torque controller over a given desired trajectory.\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.TestIfSE3-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.TestIfSE3",
    "category": "method",
    "text": "TestIfSE3(mat)\n\nReturns true if mat is close to or on the manifold SE(3).\n\nExamples\n\njulia> TestIfSE3([1.0 0.0 0.0 1.2; 0.0 0.1 -0.95 1.5; 0.0 1.0 0.1 -0.9; 0.0 0.0 0.1 0.98])\nfalse\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.TestIfSO3-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.TestIfSO3",
    "category": "method",
    "text": "TestIfSO3(mat)\n\nReturns true if mat is close to or on the manifold SO(3).\n\nExamples\n\njulia> TestIfSO3([1.0 0.0 0.0; 0.0 0.1 -0.95; 0.0 1.0 0.1])\nfalse\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.TransInv-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.TransInv",
    "category": "method",
    "text": "TransInv(T)\n\nInverts a homogeneous transformation matrix.\n\nExamples\n\njulia> TransInv([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])\n4×4 Array{Int64,2}:\n 1   0  0   0\n 0   0  1  -3\n 0  -1  0   0\n 0   0  0   1\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.TransToRp-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.TransToRp",
    "category": "method",
    "text": "TransToRp(T)\n\nConverts a homogeneous transformation matrix into a rotation matrix and position vector.\n\nExamples\n\njulia> TransToRp([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])\n([1 0 0; 0 0 -1; 0 1 0], [0, 0, 3])\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.VecTose3-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.VecTose3",
    "category": "method",
    "text": "VecTose3(V)\n\nConverts a spatial velocity vector into a 4x4 matrix in se3.\n\nExamples\n\njulia> VecTose3([1 2 3 4 5 6])\n4×4 Array{Float64,2}:\n  0.0  -3.0   2.0  4.0\n  3.0   0.0  -1.0  5.0\n -2.0   1.0   0.0  6.0\n  0.0   0.0   0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.VecToso3-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.VecToso3",
    "category": "method",
    "text": "VecToso3(omg)\n\nConverts a 3-vector to an so(3) representation.\n\nExamples\n\njulia> VecToso3([1 2 3])\n3×3 Array{Int64,2}:\n  0  -3   2\n  3   0  -1\n -2   1   0\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.VelQuadraticForces-Tuple{Array,Array,Array,Array,AbstractArray{T,2} where T}",
    "page": "Home",
    "title": "ModernRobotics.VelQuadraticForces",
    "category": "method",
    "text": "VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist)\n\nComputes the Coriolis and centripetal terms in the inverse dynamics of an open chain robot.\n\nExamples\n\njulia> VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist)\n3-element Array{Float64,1}:\n  0.26453118054501235 \n -0.0550515682891655  \n -0.006891320068248911\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.ad-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.ad",
    "category": "method",
    "text": "ad(V)\n\nCalculate the 6x6 matrix [adV] of the given 6-vector.\n\nExamples\n\njulia> ad([1, 2, 3, 4, 5, 6])\n6×6 Array{Float64,2}:\n  0.0  -3.0   2.0   0.0   0.0   0.0\n  3.0   0.0  -1.0   0.0   0.0   0.0\n -2.0   1.0   0.0   0.0   0.0   0.0\n  0.0  -6.0   5.0   0.0  -3.0   2.0\n  6.0   0.0  -4.0   3.0   0.0  -1.0\n -5.0   4.0   0.0  -2.0   1.0   0.0\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.se3ToVec-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.se3ToVec",
    "category": "method",
    "text": "se3ToVec(se3mat)\n\nConverts an se3 matrix into a spatial velocity vector.\n\nExamples\n\njulia> se3ToVec([0 -3 2 4; 3 0 -1 5; -2 1 0 6; 0 0 0 0])\n1×6 Array{Int64,2}:\n 1  2  3  4  5  6\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.so3ToVec-Tuple{Array}",
    "page": "Home",
    "title": "ModernRobotics.so3ToVec",
    "category": "method",
    "text": "so3ToVec(so3mat)\n\nConverts an so(3) representation to a 3-vector.\n\nExamples\n\njulia> so3ToVec([0 -3 2; 3 0 -1; -2 1 0])\n3-element Array{Int64,1}:\n 1\n 2\n 3\n\n\n\n\n\n"
},

{
    "location": "#ModernRobotics.jl-1",
    "page": "Home",
    "title": "ModernRobotics.jl",
    "category": "section",
    "text": "Modules = [ModernRobotics]"
},

]}

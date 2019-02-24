var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#ModernRobotics.jl-1",
    "page": "Home",
    "title": "ModernRobotics.jl",
    "category": "section",
    "text": "Some examples can be found on the Examples page.See the Index for the complete list of documented functions and types."
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The latest release of ModernRobotics can be installed from the Julia REPL prompt withjulia> Pkg.add(\"ModernRobotics\")warning: Warning\nModernRobotics.jl is currently in the process of being added to the official registry of Julia packages: JuliaLang/METADATA.jl/pull/21669.Should the above command fail please try enter Julia\'s Pkg REPL mode by pressing ] and use this one instead:(v1.1) pkg> add https://github.com/ferrolho/ModernRobotics.jl"
},

{
    "location": "#Manual-Outline-1",
    "page": "Home",
    "title": "Manual Outline",
    "category": "section",
    "text": "Pages = [\n    \"man/examples.md\",\n]\nDepth = 1"
},

{
    "location": "#Library-Outline-1",
    "page": "Home",
    "title": "Library Outline",
    "category": "section",
    "text": "Pages = [\"lib/public.md\"]"
},

{
    "location": "#main-index-1",
    "page": "Home",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"lib/public.md\"]"
},

{
    "location": "man/examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "man/examples/#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": "Here are some examples for Forward and Inverse Dynamics Trajectories."
},

{
    "location": "man/examples/#Contents-1",
    "page": "Examples",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"examples.md\"]"
},

{
    "location": "man/examples/#Inverse-Dynamics-Trajectory-1",
    "page": "Examples",
    "title": "Inverse Dynamics Trajectory",
    "category": "section",
    "text": "using ModernRobotics\n\nimport LinearAlgebra\nconst linalg = LinearAlgebra;Create a trajectory to follow using functions from Chapter 9:thetastart = [0, 0, 0]\n\nthetaend = [π / 2, π / 2, π / 2]\n\nTf = 3\n\nN = 1000\n\nmethod = 5\n\ntraj = JointTrajectory(thetastart, thetaend, Tf, N, method)\n\nthetamat = copy(traj)\ndthetamat = zeros(1000, 3)\nddthetamat = zeros(1000, 3)\n\ndt = Tf / (N - 1.0)\n\nfor i = 1:size(traj, 1) - 1\n    dthetamat[i + 1, :] = (thetamat[i + 1, :] - thetamat[i, :]) / dt\n    ddthetamat[i + 1, :] = (dthetamat[i + 1, :] - dthetamat[i, :]) / dt\nendInitialize robot description (example with 3 links):g = [0, 0, -9.8]\n\nFtipmat = ones(N, 6)\n\nM01 = [ 1 0 0        0 ;\n        0 1 0        0 ;\n        0 0 1 0.089159 ;\n        0 0 0        1 ]\n\nM12 = [ 0 0 1    0.28 ;\n        0 1 0 0.13585 ;\n       -1 0 0       0 ;\n        0 0 0       1 ]\n\nM23 = [ 1 0 0       0 ;\n        0 1 0 -0.1197 ;\n        0 0 1   0.395 ;\n        0 0 0       1 ]\n\nM34 = [ 1 0 0       0 ;\n        0 1 0       0 ;\n        0 0 1 0.14225 ;\n        0 0 0       1 ]\n\nMlist = [M01, M12, M23, M34]\n\nG1 = linalg.Diagonal([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])\nG2 = linalg.Diagonal([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])\nG3 = linalg.Diagonal([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])\n\nGlist = [G1, G2, G3]\n\nSlist = [ 1 0 1      0 1     0 ;\n          0 1 0 -0.089 0     0 ;\n          0 1 0 -0.089 0 0.425 ]\'\n\ntaumat = InverseDynamicsTrajectory(thetamat, dthetamat, ddthetamat, g, Ftipmat, Mlist, Glist, Slist)Plot the joint forces/torques:using Plots\ngr()\n\ntimestamp = range(1, Tf, length=N)\n\nplot(timestamp, taumat[:, 1], linewidth=2, label=\"Tau 1\")\nplot!(timestamp, taumat[:, 2], linewidth=2, label=\"Tau 2\")\nplot!(timestamp, taumat[:, 3], linewidth=2, label=\"Tau 3\")\n\nxlabel!(\"Time\")\nylabel!(\"Torque\")\ntitle!(\"Plot of Torque Trajectories\")(Image: inverse_dynamics_trajectory)"
},

{
    "location": "man/examples/#Forward-Dynamics-Trajectory-1",
    "page": "Examples",
    "title": "Forward Dynamics Trajectory",
    "category": "section",
    "text": "dt = 0.1\n\nintRes = 8\n\nthetalist = [0.1, 0.1, 0.1]\n\ndthetalist = [0.1, 0.2, 0.3]\n\ntaumat = [[3.63, -6.58, -5.57], [3.74, -5.55,  -5.5],\n          [4.31, -0.68, -5.19], [5.18,  5.63, -4.31],\n          [5.85,  8.17, -2.59], [5.78,  2.79,  -1.7],\n          [4.99,  -5.3, -1.19], [4.08, -9.41,  0.07],\n          [3.56, -10.1,  0.97], [3.49, -9.41,  1.23]]\n\ntaumat = cat(taumat..., dims=2)\'\n\nthetamat, dthetamat = ForwardDynamicsTrajectory(thetalist, dthetalist, taumat, g,\n                                                Ftipmat, Mlist, Glist, Slist, dt, intRes)Plot the joint angle/velocities:theta1 = thetamat[:, 1]\ntheta2 = thetamat[:, 2]\ntheta3 = thetamat[:, 3]\n\ndtheta1 = dthetamat[:, 1]\ndtheta2 = dthetamat[:, 2]\ndtheta3 = dthetamat[:, 3]\n\nN = size(taumat, 1)\nTf = size(taumat, 1) * dt\n\ntimestamp = range(0, Tf, length=N)\n\nplot(timestamp, theta1, linewidth=2, label=\"Theta1\")\nplot!(timestamp, theta2, linewidth=2, label=\"Theta2\")\nplot!(timestamp, theta3, linewidth=2, label=\"Theta3\")\n\nplot!(timestamp, dtheta1, linewidth=2, label=\"DTheta1\")\nplot!(timestamp, dtheta2, linewidth=2, label=\"DTheta2\")\nplot!(timestamp, dtheta3, linewidth=2, label=\"DTheta3\")\n\nxlabel!(\"Time\")\nylabel!(\"Joint Angles/Velocities\")\ntitle!(\"Plot of Joint Angles and Joint Velocities\")(Image: forward_dynamics_trajectory)"
},

{
    "location": "man/examples/#Simulate-Control-1",
    "page": "Examples",
    "title": "Simulate Control",
    "category": "section",
    "text": "Create a trajectory to follow:thetaend = [π / 2, π, 1.5 * π]\n\nTf = 1\ndt = 0.01\nN = round(Int, Tf / dt)\nmethod = 5\n\ntraj = JointTrajectory(thetalist, thetaend, Tf, N, method)\n\nthetamatd = copy(traj)\ndthetamatd = zeros(N, 3)\nddthetamatd = zeros(N, 3)\n\ndt = Tf / (N - 1)\n\nfor i = 1:size(traj, 1)-1\n    dthetamatd[i + 1, :] = (thetamatd[i + 1, :] - thetamatd[i, :]) / dt\n    ddthetamatd[i + 1, :] = (dthetamatd[i + 1, :] - dthetamatd[i, :]) / dt\nendCreate a (possibly) wrong robot description:gtilde = [0.8, 0.2, -8.8]\n\nMhat01 = [1 0 0   0 ;\n          0 1 0   0 ;\n          0 0 1 0.1 ;\n          0 0 0   1 ]\n\nMhat12 = [ 0 0 1 0.3 ;\n           0 1 0 0.2 ;\n          -1 0 0   0 ;\n           0 0 0   1 ]\n\nMhat23 = [1 0 0    0 ;\n          0 1 0 -0.2 ;\n          0 0 1  0.4 ;\n          0 0 0    1 ]\n\nMhat34 = [1 0 0   0 ;\n          0 1 0   0 ;\n          0 0 1 0.2 ;\n          0 0 0   1 ]\n\nMtildelist = [Mhat01, Mhat12, Mhat23, Mhat34]\n\nGhat1 = linalg.Diagonal([0.1, 0.1, 0.1, 4, 4, 4])\nGhat2 = linalg.Diagonal([0.3, 0.3, 0.1, 9, 9, 9])\nGhat3 = linalg.Diagonal([0.1, 0.1, 0.1, 3, 3, 3])\n\nGtildelist = [Ghat1, Ghat2, Ghat3]\n\nFtipmat = ones(size(traj, 1), 6)\n\nKp = 20\nKi = 10\nKd = 18\n\nintRes = 8\n\ntaumat, thetamat = SimulateControl(thetalist, dthetalist, g, Ftipmat, Mlist, Glist,\n                                   Slist, thetamatd, dthetamatd, ddthetamatd, gtilde,\n                                   Mtildelist, Gtildelist, Kp, Ki, Kd, dt, intRes)Finally, plot the results:N, links = size(thetamat)\n\nTf = N * dt\n\ntimestamp = range(0, Tf, length=N)\n\nplot()\n\nfor i = 1:links\n    plot!(timestamp, thetamat[:, i], lw=2, linestyle=:dash, label=\"ActualTheta $i\")\n    plot!(timestamp, thetamatd[:, i], lw=2, linestyle=:dot, label=\"DesiredTheta $i\")\nend\n\nxlabel!(\"Time\")\nylabel!(\"Joint Angles\")\ntitle!(\"Plot of Actual and Desired Joint Angles\")(Image: simulated_control)"
},

{
    "location": "lib/public/#",
    "page": "Public",
    "title": "Public",
    "category": "page",
    "text": ""
},

{
    "location": "lib/public/#Public-Documentation-1",
    "page": "Public",
    "title": "Public Documentation",
    "category": "section",
    "text": "Documentation for ModernRobotics.jl\'s public interface."
},

{
    "location": "lib/public/#Contents-1",
    "page": "Public",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"public.md\"]"
},

{
    "location": "lib/public/#Index-1",
    "page": "Public",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"public.md\"]"
},

{
    "location": "lib/public/#ModernRobotics.Adjoint-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.Adjoint",
    "category": "method",
    "text": "Adjoint(T)\n\nComputes the adjoint representation of a homogeneous transformation matrix.\n\nExamples\n\njulia> Adjoint([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])\n6×6 Array{Float64,2}:\n 1.0  0.0   0.0  0.0  0.0   0.0\n 0.0  0.0  -1.0  0.0  0.0   0.0\n 0.0  1.0   0.0  0.0  0.0   0.0\n 0.0  0.0   3.0  1.0  0.0   0.0\n 3.0  0.0   0.0  0.0  0.0  -1.0\n 0.0  0.0   0.0  0.0  1.0   0.0\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.AxisAng3-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.AxisAng3",
    "category": "method",
    "text": "AxisAng3(expc3)\n\nConverts a 3-vector of exponential coordinates for rotation into axis-angle form.\n\nExamples\n\njulia> AxisAng3([1, 2, 3])\n([0.267261, 0.534522, 0.801784], 3.7416573867739413)\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.AxisAng6-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.AxisAng6",
    "category": "method",
    "text": "AxisAng6(expc6)\n\nConverts a 6-vector of exponential coordinates into screw axis-angle form.\n\nExamples\n\njulia> AxisAng6([1, 0, 0, 1, 2, 3])\n([1.0, 0.0, 0.0, 1.0, 2.0, 3.0], 1.0)\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.ComputedTorque-Tuple{Array,Array,Array,Array,Array,Array,AbstractArray{T,2} where T,Array,Array,Array,Number,Number,Number}",
    "page": "Public",
    "title": "ModernRobotics.ComputedTorque",
    "category": "method",
    "text": "ComputedTorque(thetalist, dthetalist, eint, g, Mlist, Glist, Slist, thetalistd, dthetalistd, ddthetalistd, Kp, Ki, Kd)\n\nComputes the joint control torques at a particular time instant.\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.CubicTimeScaling-Tuple{Number,Number}",
    "page": "Public",
    "title": "ModernRobotics.CubicTimeScaling",
    "category": "method",
    "text": "CubicTimeScaling(Tf, t)\n\nComputes s(t) for a cubic time scaling.\n\nExamples\n\njulia> CubicTimeScaling(2, 0.6)\n0.21600000000000003\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.DistanceToSE3-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.DistanceToSE3",
    "category": "method",
    "text": "DistanceToSE3(mat)\n\nReturns the Frobenius norm to describe the distance of mat from the SE(3) manifold.\n\nExamples\n\njulia> DistanceToSE3([1.0 0.0 0.0 1.2; 0.0 0.1 -0.95 1.5; 0.0 1.0 0.1 -0.9; 0.0 0.0 0.1 0.98])\n0.13493053768513638\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.DistanceToSO3-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.DistanceToSO3",
    "category": "method",
    "text": "DistanceToSO3(mat)\n\nReturns the Frobenius norm to describe the distance of mat from the SO(3) manifold.\n\nExamples\n\njulia> DistanceToSO3([1.0 0.0 0.0; 0.0 0.1 -0.95; 0.0 1.0 0.1])\n0.08835298523536149\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.EndEffectorForces-Tuple{Array,Array,Array,Array,AbstractArray{T,2} where T}",
    "page": "Public",
    "title": "ModernRobotics.EndEffectorForces",
    "category": "method",
    "text": "EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist)\n\nComputes the joint forces/torques an open chain robot requires only to create the end-effector force Ftip.\n\nArguments\n\nthetalist: the n-vector of joint variables.\nFtip: the spatial force applied by the end-effector expressed in frame {n+1}.\nMlist: the list of link frames i relative to i-1 at the home position.\nGlist: the spatial inertia matrices Gi of the links.\nSlist: the screw axes Si of the joints in a space frame, in the format of a matrix with axes as the columns.\n\nReturns the joint forces and torques required only to create the end-effector force Ftip. This function calls InverseDynamics with g = 0, dthetalist = 0, and ddthetalist = 0.\n\nExamples\n\njulia> import LinearAlgebra\n\njulia> const linalg = LinearAlgebra;\n\njulia> thetalist = [0.1, 0.1, 0.1]\n3-element Array{Float64,1}:\n 0.1\n 0.1\n 0.1\n\njulia> Ftip = [1, 1, 1, 1, 1, 1]\n6-element Array{Int64,1}:\n 1\n 1\n 1\n 1\n 1\n 1\n\njulia> M01 = [1 0 0        0;\n              0 1 0        0;\n              0 0 1 0.089159;\n              0 0 0        1]\n4×4 Array{Float64,2}:\n 1.0  0.0  0.0  0.0     \n 0.0  1.0  0.0  0.0     \n 0.0  0.0  1.0  0.089159\n 0.0  0.0  0.0  1.0     \n\njulia> M12 = [ 0 0 1    0.28;\n               0 1 0 0.13585;\n              -1 0 0       0;\n               0 0 0       1]\n4×4 Array{Float64,2}:\n  0.0  0.0  1.0  0.28   \n  0.0  1.0  0.0  0.13585\n -1.0  0.0  0.0  0.0    \n  0.0  0.0  0.0  1.0    \n\njulia> M23 = [1 0 0       0;\n              0 1 0 -0.1197;\n              0 0 1   0.395;\n              0 0 0       1]\n4×4 Array{Float64,2}:\n 1.0  0.0  0.0   0.0   \n 0.0  1.0  0.0  -0.1197\n 0.0  0.0  1.0   0.395 \n 0.0  0.0  0.0   1.0   \n\njulia> M34 = [1 0 0       0;\n              0 1 0       0;\n              0 0 1 0.14225;\n              0 0 0       1]\n4×4 Array{Float64,2}:\n 1.0  0.0  0.0  0.0    \n 0.0  1.0  0.0  0.0    \n 0.0  0.0  1.0  0.14225\n 0.0  0.0  0.0  1.0    \n\njulia> Mlist = [M01, M12, M23, M34]\n4-element Array{Array{Float64,2},1}:\n [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.089159; 0.0 0.0 0.0 1.0] \n [0.0 0.0 1.0 0.28; 0.0 1.0 0.0 0.13585; -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0]\n [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 -0.1197; 0.0 0.0 1.0 0.395; 0.0 0.0 0.0 1.0]\n [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.14225; 0.0 0.0 0.0 1.0]  \n\njulia> G1 = linalg.Diagonal([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])\n6×6 LinearAlgebra.Diagonal{Float64,Array{Float64,1}}:\n 0.010267   ⋅         ⋅        ⋅    ⋅    ⋅ \n  ⋅        0.010267   ⋅        ⋅    ⋅    ⋅ \n  ⋅         ⋅        0.00666   ⋅    ⋅    ⋅ \n  ⋅         ⋅         ⋅       3.7   ⋅    ⋅ \n  ⋅         ⋅         ⋅        ⋅   3.7   ⋅ \n  ⋅         ⋅         ⋅        ⋅    ⋅   3.7\n\njulia> G2 = linalg.Diagonal([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])\n6×6 LinearAlgebra.Diagonal{Float64,Array{Float64,1}}:\n 0.22689   ⋅        ⋅          ⋅      ⋅      ⋅   \n  ⋅       0.22689   ⋅          ⋅      ⋅      ⋅   \n  ⋅        ⋅       0.0151074   ⋅      ⋅      ⋅   \n  ⋅        ⋅        ⋅         8.393   ⋅      ⋅   \n  ⋅        ⋅        ⋅          ⋅     8.393   ⋅   \n  ⋅        ⋅        ⋅          ⋅      ⋅     8.393\n\njulia> G3 = linalg.Diagonal([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])\n6×6 LinearAlgebra.Diagonal{Float64,Array{Float64,1}}:\n 0.0494433   ⋅          ⋅         ⋅      ⋅      ⋅   \n  ⋅         0.0494433   ⋅         ⋅      ⋅      ⋅   \n  ⋅          ⋅         0.004095   ⋅      ⋅      ⋅   \n  ⋅          ⋅          ⋅        2.275   ⋅      ⋅   \n  ⋅          ⋅          ⋅         ⋅     2.275   ⋅   \n  ⋅          ⋅          ⋅         ⋅      ⋅     2.275\n\njulia> Glist = [G1, G2, G3]\n3-element Array{LinearAlgebra.Diagonal{Float64,Array{Float64,1}},1}:\n [0.010267 0.0 … 0.0 0.0; 0.0 0.010267 … 0.0 0.0; … ; 0.0 0.0 … 3.7 0.0; 0.0 0.0 … 0.0 3.7]      \n [0.22689 0.0 … 0.0 0.0; 0.0 0.22689 … 0.0 0.0; … ; 0.0 0.0 … 8.393 0.0; 0.0 0.0 … 0.0 8.393]    \n [0.0494433 0.0 … 0.0 0.0; 0.0 0.0494433 … 0.0 0.0; … ; 0.0 0.0 … 2.275 0.0; 0.0 0.0 … 0.0 2.275]\n\njulia> Slist = [ 1  0  1      0  1      0;\n                 0  1  0 -0.089  0      0;\n                 0  1  0 -0.089  0  0.425]\'\n6×3 LinearAlgebra.Adjoint{Float64,Array{Float64,2}}:\n 1.0   0.0     0.0  \n 0.0   1.0     1.0  \n 1.0   0.0     0.0  \n 0.0  -0.089  -0.089\n 1.0   0.0     0.0  \n 0.0   0.0     0.425\n\njulia> EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist)\n3-element Array{Float64,1}:\n 1.4095460782639782\n 1.8577149723180628\n 1.392409          \n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.EulerStep-Tuple{Array,Array,Array,Number}",
    "page": "Public",
    "title": "ModernRobotics.EulerStep",
    "category": "method",
    "text": "EulerStep(thetalist, dthetalist, ddthetalist, dt)\n\nCompute the joint angles and velocities at the next timestep using first order Euler integration.\n\nArguments\n\nthetalist: the n-vector of joint variables.\ndthetalist: the n-vector of joint rates.\nddthetalist: the n-vector of joint accelerations.\ndt: the timestep delta t.\n\nReturn\n\nthetalistNext: the vector of joint variables after dt from first order Euler integration.\ndthetalistNext: the vector of joint rates after dt from first order Euler integration.\n\nExamples\n\njulia> EulerStep([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], [2, 1.5, 1], 0.1)\n([0.11, 0.12, 0.13], [0.3, 0.35, 0.4])\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.FKinBody-Tuple{AbstractArray{T,2} where T,AbstractArray{T,2} where T,Array}",
    "page": "Public",
    "title": "ModernRobotics.FKinBody",
    "category": "method",
    "text": "FKinBody(M, Blist, thetalist)\n\nComputes forward kinematics in the body frame for an open chain robot.\n\nExamples\n\njulia> M = [ -1  0  0  0 ;\n              0  1  0  6 ;\n              0  0 -1  2 ;\n              0  0  0  1 ];\n\njulia> Blist = [  0  0 -1  2  0  0   ;\n                  0  0  0  0  1  0   ;\n                  0  0  1  0  0  0.1 ]\';\n\njulia> thetalist = [ π/2, 3, π ];\n\njulia> FKinBody(M, Blist, thetalist)\n4×4 Array{Float64,2}:\n -1.14424e-17  1.0           0.0  -5.0    \n  1.0          1.14424e-17   0.0   4.0    \n  0.0          0.0          -1.0   1.68584\n  0.0          0.0           0.0   1.0    \n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.FKinSpace-Tuple{AbstractArray{T,2} where T,AbstractArray{T,2} where T,Array}",
    "page": "Public",
    "title": "ModernRobotics.FKinSpace",
    "category": "method",
    "text": "FKinSpace(M, Slist, thetalist)\n\nComputes forward kinematics in the space frame for an open chain robot.\n\nExamples\n\njulia> M = [ -1  0  0  0 ;\n              0  1  0  6 ;\n              0  0 -1  2 ;\n              0  0  0  1 ];\n\njulia> Slist = [  0  0  1  4  0  0   ;\n                  0  0  0  0  1  0   ;\n                  0  0 -1 -6  0 -0.1 ]\';\n\njulia> thetalist = [ π/2, 3, π ];\n\njulia> FKinSpace(M, Slist, thetalist)\n4×4 Array{Float64,2}:\n -1.14424e-17  1.0           0.0  -5.0    \n  1.0          1.14424e-17   0.0   4.0    \n  0.0          0.0          -1.0   1.68584\n  0.0          0.0           0.0   1.0    \n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.ForwardDynamics-Tuple{Array,Array,Array,Array,Array,Array,Array,AbstractArray{T,2} where T}",
    "page": "Public",
    "title": "ModernRobotics.ForwardDynamics",
    "category": "method",
    "text": "ForwardDynamics(thetalist, dthetalist, taulist, g, Ftip, Mlist, Glist, Slist)\n\nComputes forward dynamics in the space frame for an open chain robot.\n\nExamples\n\njulia> ForwardDynamics(thetalist, dthetalist, taulist, g, Ftip, Mlist, Glist, Slist)\n3-element Array{Float64,1}:\n  -0.9739290670855626\n  25.584667840340558 \n -32.91499212478149  \n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.ForwardDynamicsTrajectory-Tuple{Array,Array,AbstractArray{T,2} where T,Array,Array,Array,Array,AbstractArray{T,2} where T,Number,Number}",
    "page": "Public",
    "title": "ModernRobotics.ForwardDynamicsTrajectory",
    "category": "method",
    "text": "ForwardDynamicsTrajectory(thetalist, dthetalist, taumat, g, Ftipmat, Mlist, Glist, Slist, dt, intRes)\n\nSimulates the motion of a serial chain given an open-loop history of joint forces/torques.\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.GravityForces-Tuple{Array,Array,Array,Array,AbstractArray{T,2} where T}",
    "page": "Public",
    "title": "ModernRobotics.GravityForces",
    "category": "method",
    "text": "GravityForces(thetalist, g, Mlist, Glist, Slist)\n\nComputes the joint forces/torques an open chain robot requires to overcome gravity at its configuration.\n\nExamples\n\njulia> GravityForces(thetalist, g, Mlist, Glist, Slist)\n3-element Array{Float64,1}:\n  28.40331261821983  \n -37.64094817177068  \n  -5.4415891999683605\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.IKinBody-Tuple{AbstractArray{T,2} where T,AbstractArray{T,2} where T,AbstractArray{T,2} where T,Array,Number,Number}",
    "page": "Public",
    "title": "ModernRobotics.IKinBody",
    "category": "method",
    "text": "IKinBody(Blist, M, T, thetalist0, eomg, ev)\n\nComputes inverse kinematics in the body frame for an open chain robot.\n\nExamples\n\njulia> Blist = [  0  0 -1  2  0  0   ;\n                  0  0  0  0  1  0   ;\n                  0  0  1  0  0  0.1 ]\';\n\njulia> M = [ -1  0  0  0 ;\n              0  1  0  6 ;\n              0  0 -1  2 ;\n              0  0  0  1 ];\n\njulia> T = [ 0  1  0     -5 ;\n             1  0  0      4 ;\n             0  0 -1 1.6858 ;\n             0  0  0      1 ];\n\njulia> thetalist0 = [1.5, 2.5, 3];\n\njulia> eomg, ev = 0.01, 0.001;\n\njulia> IKinBody(Blist, M, T, thetalist0, eomg, ev)\n([1.57074, 2.99967, 3.14154], true)\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.IKinSpace-Tuple{AbstractArray{T,2} where T,AbstractArray{T,2} where T,AbstractArray{T,2} where T,Array,Number,Number}",
    "page": "Public",
    "title": "ModernRobotics.IKinSpace",
    "category": "method",
    "text": "IKinSpace(Slist, M, T, thetalist0, eomg, ev)\n\nComputes inverse kinematics in the space frame for an open chain robot.\n\nExamples\n\njulia> Slist = [  0  0  1  4  0  0   ;\n                  0  0  0  0  1  0   ;\n                  0  0 -1 -6  0 -0.1 ]\';\n\njulia> M = [ -1  0  0  0 ;\n              0  1  0  6 ;\n              0  0 -1  2 ;\n              0  0  0  1 ];\n\njulia> T = [ 0  1  0     -5 ;\n             1  0  0      4 ;\n             0  0 -1 1.6858 ;\n             0  0  0      1 ];\n\njulia> thetalist0 = [1.5, 2.5, 3];\n\njulia> eomg, ev = 0.01, 0.001;\n\njulia> IKinSpace(Slist, M, T, thetalist0, eomg, ev)\n([1.57074, 2.99966, 3.14153], true)\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.InverseDynamics-Tuple{Array,Array,Array,Array,Array,Array,Array,AbstractArray{T,2} where T}",
    "page": "Public",
    "title": "ModernRobotics.InverseDynamics",
    "category": "method",
    "text": "InverseDynamics(thetalist, dthetalist, ddthetalist, g, Ftip, Mlist, Glist, Slist)\n\nComputes inverse dynamics in the space frame for an open chain robot.\n\nExamples\n\njulia> InverseDynamics(thetalist, dthetalist, ddthetalist, g, Ftip, Mlist, Glist, Slist)\n3-element Array{Float64,1}:\n  74.69616155287451 \n -33.06766015851458 \n  -3.230573137901424\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.InverseDynamicsTrajectory-Tuple{Array,Array,Array,Array,Array,Array,Array,AbstractArray{T,2} where T}",
    "page": "Public",
    "title": "ModernRobotics.InverseDynamicsTrajectory",
    "category": "method",
    "text": "InverseDynamicsTrajectory(thetamat, dthetamat, ddthetamat, g, Ftipmat, Mlist, Glist, Slist)\n\nCalculates the joint forces/torques required to move the serial chain along the given trajectory using inverse dynamics.\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.JacobianBody-Tuple{AbstractArray{T,2} where T,Array}",
    "page": "Public",
    "title": "ModernRobotics.JacobianBody",
    "category": "method",
    "text": "JacobianBody(Blist, thetalist)\n\nComputes the body Jacobian for an open chain robot.\n\nExamples\n\njulia> Blist = [0 0 1   0 0.2 0.2;\n                1 0 0   2   0   3;\n                0 1 0   0   2   1;\n                1 0 0 0.2 0.3 0.4]\';\n\njulia> thetalist = [0.2, 1.1, 0.1, 1.2];\n\njulia> JacobianBody(Blist, thetalist)\n6×4 Array{Float64,2}:\n -0.0452841  0.995004    0.0       1.0\n  0.743593   0.0930486   0.362358  0.0\n -0.667097   0.0361754  -0.932039  0.0\n  2.32586    1.66809     0.564108  0.2\n -1.44321    2.94561     1.43307   0.3\n -2.0664     1.82882    -1.58869   0.4\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.JacobianSpace-Tuple{AbstractArray{T,2} where T,Array}",
    "page": "Public",
    "title": "ModernRobotics.JacobianSpace",
    "category": "method",
    "text": "JacobianSpace(Slist, thetalist)\n\nComputes the space Jacobian for an open chain robot.\n\nExamples\n\njulia> Slist = [0 0 1   0 0.2 0.2;\n                1 0 0   2   0   3;\n                0 1 0   0   2   1;\n                1 0 0 0.2 0.3 0.4]\';\n\njulia> thetalist = [0.2, 1.1, 0.1, 1.2];\n\njulia> JacobianSpace(Slist, thetalist)\n6×4 Array{Float64,2}:\n 0.0  0.980067  -0.0901156   0.957494 \n 0.0  0.198669   0.444554    0.284876 \n 1.0  0.0        0.891207   -0.0452841\n 0.0  1.95219   -2.21635    -0.511615 \n 0.2  0.436541  -2.43713     2.77536  \n 0.2  2.96027    3.23573     2.22512  \n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.JointTrajectory-Tuple{Array,Array,Number,Integer,Integer}",
    "page": "Public",
    "title": "ModernRobotics.JointTrajectory",
    "category": "method",
    "text": "JointTrajectory(thetastart, thetaend, Tf, N, method)\n\nComputes a straight-line trajectory in joint space.\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.MassMatrix-Tuple{Array,Array,Array,AbstractArray{T,2} where T}",
    "page": "Public",
    "title": "ModernRobotics.MassMatrix",
    "category": "method",
    "text": "MassMatrix(thetalist, Mlist, Glist, Slist)\n\nComputes the mass matrix of an open chain robot based on the given configuration.\n\nExamples\n\njulia> MassMatrix(thetalist, Mlist, Glist, Slist)\n3×3 Array{Float64,2}:\n 22.5433      -0.307147  -0.00718426\n -0.307147     1.96851    0.432157  \n -0.00718426   0.432157   0.191631  \n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.MatrixExp3-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.MatrixExp3",
    "category": "method",
    "text": "MatrixExp3(so3mat)\n\nComputes the matrix exponential of a matrix in so(3).\n\nExamples\n\njulia> MatrixExp3([0 -3 2; 3 0 -1; -2 1 0])\n3×3 Array{Float64,2}:\n -0.694921   0.713521  0.0892929\n -0.192007  -0.303785  0.933192 \n  0.692978   0.63135   0.348107 \n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.MatrixExp6-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.MatrixExp6",
    "category": "method",
    "text": "MatrixExp6(se3mat)\n\nComputes the matrix exponential of an se3 representation of exponential coordinates.\n\nExamples\n\njulia> MatrixExp6([0 0 0 0; 0 0 -1.57079632 2.35619449; 0 1.57079632 0 2.35619449; 0 0 0 0])\n4×4 Array{Float64,2}:\n 1.0  0.0         0.0        0.0       \n 0.0  6.7949e-9  -1.0        1.01923e-8\n 0.0  1.0         6.7949e-9  3.0       \n 0.0  0.0         0.0        1.0       \n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.MatrixLog3-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.MatrixLog3",
    "category": "method",
    "text": "MatrixLog3(R)\n\nComputes the matrix logarithm of a rotation matrix.\n\nExamples\n\njulia> MatrixLog3([0 0 1; 1 0 0; 0 1 0])\n3×3 Array{Float64,2}:\n  0.0     -1.2092   1.2092\n  1.2092   0.0     -1.2092\n -1.2092   1.2092   0.0   \n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.MatrixLog6-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.MatrixLog6",
    "category": "method",
    "text": "MatrixLog6(T)\n\nComputes the matrix logarithm of a homogeneous transformation matrix.\n\nExamples\n\njulia> MatrixLog6([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])\n4×4 Array{Float64,2}:\n 0.0  0.0      0.0     0.0    \n 0.0  0.0     -1.5708  2.35619\n 0.0  1.5708   0.0     2.35619\n 0.0  0.0      0.0     0.0    \n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.NearZero-Tuple{Number}",
    "page": "Public",
    "title": "ModernRobotics.NearZero",
    "category": "method",
    "text": "NearZero(z)\n\nDetermines whether a scalar is small enough to be treated as zero.\n\nExamples\n\njulia> NearZero(-1e-7)\ntrue\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.Normalize-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.Normalize",
    "category": "method",
    "text": "Normalize(V)\n\nNormalizes a vector.\n\nExamples\n\njulia> Normalize([1, 2, 3])\n3-element Array{Float64,1}:\n 0.2672612419124244\n 0.5345224838248488\n 0.8017837257372732\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.ProjectToSE3-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.ProjectToSE3",
    "category": "method",
    "text": "ProjectToSE3(mat)\n\nReturns a projection of mat into SE(3).\n\nExamples\n\njulia> ProjectToSE3([0.675 0.150 0.720 1.2; 0.370 0.771 -0.511 5.4; -0.630 0.619 0.472 3.6; 0.003 0.002 0.010 0.9])\n4×4 Array{Float64,2}:\n  0.679011  0.148945   0.718859  1.2\n  0.373207  0.773196  -0.512723  5.4\n -0.632187  0.616428   0.469421  3.6\n  0.0       0.0        0.0       1.0\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.ProjectToSO3-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.ProjectToSO3",
    "category": "method",
    "text": "ProjectToSO3(mat)\n\nReturns a projection of mat into SO(3).\n\nExamples\n\njulia> ProjectToSO3([0.675 0.150  0.720; 0.370 0.771 -0.511; -0.630 0.619  0.472])\n3×3 Array{Float64,2}:\n  0.679011  0.148945   0.718859\n  0.373207  0.773196  -0.512723\n -0.632187  0.616428   0.469421\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.QuinticTimeScaling-Tuple{Number,Number}",
    "page": "Public",
    "title": "ModernRobotics.QuinticTimeScaling",
    "category": "method",
    "text": "QuinticTimeScaling(Tf, t)\n\nComputes s(t) for a quintic time scaling.\n\nExamples\n\njulia> QuinticTimeScaling(2, 0.6)\n0.16308\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.RotInv-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.RotInv",
    "category": "method",
    "text": "RotInv(R)\n\nInverts a rotation matrix.\n\nExamples\n\njulia> RotInv([0 0 1; 1 0 0; 0 1 0])\n3×3 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:\n 0  1  0\n 0  0  1\n 1  0  0\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.RpToTrans-Tuple{Array,Array}",
    "page": "Public",
    "title": "ModernRobotics.RpToTrans",
    "category": "method",
    "text": "RpToTrans(R, p)\n\nConverts a rotation matrix and a position vector into homogeneous transformation matrix.\n\nExamples\n\njulia> RpToTrans([1 0 0; 0 0 -1; 0 1 0], [1, 2, 5])\n4×4 Array{Int64,2}:\n 1  0   0  1\n 0  0  -1  2\n 0  1   0  5\n 0  0   0  1\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.ScrewToAxis-Tuple{Array,Array,Number}",
    "page": "Public",
    "title": "ModernRobotics.ScrewToAxis",
    "category": "method",
    "text": "ScrewToAxis(q, s, h)\n\nTakes a parametric description of a screw axis and converts it to a normalized screw axis.\n\nExamples\n\njulia> ScrewToAxis([3; 0; 0], [0; 0; 1], 2)\n6-element Array{Int64,1}:\n  0\n  0\n  1\n  0\n -3\n  2\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.SimulateControl-Tuple{Array,Array,Array,Array,Array,Array,AbstractArray{T,2} where T,Array,Array,Array,Array,Array,Array,Number,Number,Number,Number,Number}",
    "page": "Public",
    "title": "ModernRobotics.SimulateControl",
    "category": "method",
    "text": "SimulateControl(thetalist, dthetalist, g, Ftipmat, Mlist, Glist,\n                Slist, thetamatd, dthetamatd, ddthetamatd, gtilde,\n                Mtildelist, Gtildelist, Kp, Ki, Kd, dt, intRes)\n\nSimulates the computed torque controller over a given desired trajectory.\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.TestIfSE3-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.TestIfSE3",
    "category": "method",
    "text": "TestIfSE3(mat)\n\nReturns true if mat is close to or on the manifold SE(3).\n\nExamples\n\njulia> TestIfSE3([1.0 0.0 0.0 1.2; 0.0 0.1 -0.95 1.5; 0.0 1.0 0.1 -0.9; 0.0 0.0 0.1 0.98])\nfalse\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.TestIfSO3-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.TestIfSO3",
    "category": "method",
    "text": "TestIfSO3(mat)\n\nReturns true if mat is close to or on the manifold SO(3).\n\nExamples\n\njulia> TestIfSO3([1.0 0.0 0.0; 0.0 0.1 -0.95; 0.0 1.0 0.1])\nfalse\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.TransInv-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.TransInv",
    "category": "method",
    "text": "TransInv(T)\n\nInverts a homogeneous transformation matrix.\n\nExamples\n\njulia> TransInv([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])\n4×4 Array{Int64,2}:\n 1   0  0   0\n 0   0  1  -3\n 0  -1  0   0\n 0   0  0   1\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.TransToRp-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.TransToRp",
    "category": "method",
    "text": "TransToRp(T)\n\nConverts a homogeneous transformation matrix into a rotation matrix and position vector.\n\nExamples\n\njulia> TransToRp([1 0 0 0; 0 0 -1 0; 0 1 0 3; 0 0 0 1])\n([1 0 0; 0 0 -1; 0 1 0], [0, 0, 3])\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.VecTose3-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.VecTose3",
    "category": "method",
    "text": "VecTose3(V)\n\nConverts a spatial velocity vector into a 4x4 matrix in se3.\n\nExamples\n\njulia> VecTose3([1 2 3 4 5 6])\n4×4 Array{Float64,2}:\n  0.0  -3.0   2.0  4.0\n  3.0   0.0  -1.0  5.0\n -2.0   1.0   0.0  6.0\n  0.0   0.0   0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.VecToso3-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.VecToso3",
    "category": "method",
    "text": "VecToso3(omg)\n\nConverts a 3-vector to an so(3) representation.\n\nExamples\n\njulia> VecToso3([1 2 3])\n3×3 Array{Int64,2}:\n  0  -3   2\n  3   0  -1\n -2   1   0\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.VelQuadraticForces-Tuple{Array,Array,Array,Array,AbstractArray{T,2} where T}",
    "page": "Public",
    "title": "ModernRobotics.VelQuadraticForces",
    "category": "method",
    "text": "VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist)\n\nComputes the Coriolis and centripetal terms in the inverse dynamics of an open chain robot.\n\nExamples\n\njulia> VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist)\n3-element Array{Float64,1}:\n  0.26453118054501235 \n -0.0550515682891655  \n -0.006891320068248911\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.ad-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.ad",
    "category": "method",
    "text": "ad(V)\n\nCalculate the 6x6 matrix [adV] of the given 6-vector.\n\nExamples\n\njulia> ad([1, 2, 3, 4, 5, 6])\n6×6 Array{Float64,2}:\n  0.0  -3.0   2.0   0.0   0.0   0.0\n  3.0   0.0  -1.0   0.0   0.0   0.0\n -2.0   1.0   0.0   0.0   0.0   0.0\n  0.0  -6.0   5.0   0.0  -3.0   2.0\n  6.0   0.0  -4.0   3.0   0.0  -1.0\n -5.0   4.0   0.0  -2.0   1.0   0.0\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.se3ToVec-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.se3ToVec",
    "category": "method",
    "text": "se3ToVec(se3mat)\n\nConverts an se3 matrix into a spatial velocity vector.\n\nExamples\n\njulia> se3ToVec([0 -3 2 4; 3 0 -1 5; -2 1 0 6; 0 0 0 0])\n6-element Array{Int64,1}:\n 1\n 2\n 3\n 4\n 5\n 6\n\n\n\n\n\n"
},

{
    "location": "lib/public/#ModernRobotics.so3ToVec-Tuple{Array}",
    "page": "Public",
    "title": "ModernRobotics.so3ToVec",
    "category": "method",
    "text": "so3ToVec(so3mat)\n\nConverts an so(3) representation to a 3-vector.\n\nExamples\n\njulia> so3ToVec([0 -3 2; 3 0 -1; -2 1 0])\n3-element Array{Int64,1}:\n 1\n 2\n 3\n\n\n\n\n\n"
},

{
    "location": "lib/public/#Public-Interface-1",
    "page": "Public",
    "title": "Public Interface",
    "category": "section",
    "text": "Modules = [ModernRobotics]"
},

]}

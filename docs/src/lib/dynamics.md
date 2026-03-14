# Dynamics of Open Chains

*Textbook Chapter 8*

```@autodocs
Modules = [ModernRoboticsBook]
Pages = ["ModernRoboticsBook.jl"]
Filter = t -> nameof(t) in [
    :ad, :InverseDynamics, :MassMatrix, :VelQuadraticForces,
    :GravityForces, :EndEffectorForces, :ForwardDynamics,
    :EulerStep, :InverseDynamicsTrajectory, :ForwardDynamicsTrajectory,
]
```

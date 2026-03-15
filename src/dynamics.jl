export ad,
    inverse_dynamics_rnea,
    inverse_dynamics_rnea!,
    mass_matrix_crba,
    mass_matrix_crba!,
    mass_matrix_rnea,
    velocity_quadratic_forces,
    velocity_quadratic_forces!,
    gravity_forces,
    gravity_forces!,
    end_effector_forces,
    end_effector_forces!,
    forward_dynamics_aba,
    forward_dynamics_aba!,
    forward_dynamics_crba,
    forward_dynamics_crba!,
    forward_dynamics_rnea,
    forward_dynamics_rnea!,
    euler_step,
    inverse_dynamics_trajectory,
    forward_dynamics_trajectory

"""
    ad(V)

Computes the ``6 \\times 6`` matrix ``[\\text{ad}_V]`` used to calculate the Lie bracket ``[V_1, V_2] = [\\text{ad}_{V_1}] V_2``.

!!! info "What is the Lie bracket?"
    The Lie bracket measures how two twists interact — it captures the velocity
    produced by the non-commutativity of two motions. In dynamics, the ``[\\text{ad}]``
    matrix appears in the recursive Newton-Euler algorithm, where it accounts for
    how the velocity of one link affects the forces on adjacent links through
    Coriolis and centripetal effects.

# Arguments
- `V`: a 6-vector spatial velocity (twist).

# Returns
The ``6 \\times 6`` matrix ``[\\text{ad}_V]``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> ad([1, 2, 3, 4, 5, 6])
6×6 StaticArraysCore.SMatrix{6, 6, Int64, 36} with indices SOneTo(6)×SOneTo(6):
  0  -3   2   0   0   0
  3   0  -1   0   0   0
 -2   1   0   0   0   0
  0  -6   5   0  -3   2
  6   0  -4   3   0  -1
 -5   4   0  -2   1   0
```
"""
function ad(V::AbstractVector)
    ω1, ω2, ω3, v1, v2, v3 = V
    z = zero(ω1)
    @SMatrix [
        z -ω3 ω2 z z z
        ω3 z -ω1 z z z
        -ω2 ω1 z z z z
        z -v3 v2 z -ω3 ω2
        v3 z -v1 ω3 z -ω1
        -v2 v1 z -ω2 ω1 z
    ]
end

"""
    inverse_dynamics_rnea(joint_positions, joint_velocities, joint_accelerations, gravity, tip_wrench, link_frames, spatial_inertias, screw_axes)

Computes inverse dynamics in the space frame for an open chain robot using
forward-backward Newton-Euler iterations.

!!! info "Recursive Newton-Euler algorithm"
    Inverse dynamics answers: "what joint torques do I need to produce a desired motion?"
    It uses the recursive Newton-Euler algorithm, which works in two passes:

    1. **Forward pass** (base → tip): propagate link velocities and accelerations
       outward along the chain, computing each link's twist and acceleration from its
       parent's motion plus the joint's contribution.
    2. **Backward pass** (tip → base): propagate wrenches inward, computing the force
       each link needs (from its inertia and acceleration) and projecting onto joint
       axes to get the required torque.

    This is ``O(n)`` in the number of joints — much faster than forming and inverting
    the full dynamics equation. It is also the workhorse behind most dynamics
    computations: [`mass_matrix_crba`](@ref), [`velocity_quadratic_forces`](@ref),
    [`gravity_forces`](@ref), and [`end_effector_forces`](@ref) are all computed by
    calling this function with specific inputs zeroed out.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `joint_velocities`: the ``n``-vector of joint rates.
- `joint_accelerations`: the ``n``-vector of joint accelerations.
- `gravity`: the gravity vector ``g`` (e.g., `[0, 0, -9.8]`).
- `tip_wrench`: the wrench ``\\mathcal{F}_{\\text{tip}}`` applied by the end-effector expressed in frame ``\\{n+1\\}``.
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}`` (the end-effector frame relative to the last link).
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.

# Returns
The ``n``-vector of required joint forces/torques.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> inverse_dynamics_rnea([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], [2, 1.5, 1], [0, 0, -9.8], [1, 1, 1, 1, 1, 1], link_frames, spatial_inertias, screw_axes)
3-element Vector{Float64}:
  74.6961615528745
 -33.067660158514585
  -3.2305731379014246
```
"""
function inverse_dynamics_rnea(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_accelerations::AbstractVector,
    gravity::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    joint_torques = zeros(n)
    Ai = Vector{SVector{6,Float64}}(undef, n)
    AdTi = Vector{SMatrix{6,6,Float64,36}}(undef, n + 1)
    Vi = Vector{SVector{6,Float64}}(undef, n + 1)
    Vdi = Vector{SVector{6,Float64}}(undef, n + 1)
    inverse_dynamics_rnea!(
        joint_torques,
        joint_positions,
        joint_velocities,
        joint_accelerations,
        gravity,
        tip_wrench,
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Vi,
        Vdi,
    )
end

"""
    inverse_dynamics_rnea!(joint_torques, joint_positions, joint_velocities, joint_accelerations, gravity, tip_wrench, link_frames, spatial_inertias, screw_axes, Ai, AdTi, Vi, Vdi)

In-place version of [`inverse_dynamics_rnea`](@ref). Writes the result into `joint_torques`.
The workspace vectors `Ai` (length n), `AdTi` (length n+1), `Vi` (length n+1),
and `Vdi` (length n+1) must be pre-allocated.
"""
function inverse_dynamics_rnea!(
    joint_torques::AbstractVector,
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_accelerations::AbstractVector,
    gravity::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
    Ai::AbstractVector{SVector{6,Float64}},
    AdTi::AbstractVector{SMatrix{6,6,Float64,36}},
    Vi::AbstractVector{SVector{6,Float64}},
    Vdi::AbstractVector{SVector{6,Float64}},
)
    n = length(joint_positions)
    Mi = SMatrix{4,4,Float64}(LA.I)
    Vi[1] = @SVector zeros(6)
    Vdi[1] = SA[0.0, 0.0, 0.0, -gravity[1], -gravity[2], -gravity[3]]
    AdTi[n+1] = adjoint_representation(transform_inv(link_frames[n+1]))

    for i = 1:n
        Mi = Mi * SMatrix{4,4}(link_frames[i])
        Ai[i] =
            adjoint_representation(transform_inv(Mi)) * SVector{6}(@view screw_axes[:, i])
        AdTi[i] = adjoint_representation(
            matrix_exp6(vec_to_se3(Ai[i] * -joint_positions[i])) *
            transform_inv(SMatrix{4,4}(link_frames[i])),
        )
        Vi[i+1] = AdTi[i] * Vi[i] + Ai[i] * joint_velocities[i]
        Vdi[i+1] =
            AdTi[i] * Vdi[i] +
            Ai[i] * joint_accelerations[i] +
            ad(Vi[i+1]) * Ai[i] * joint_velocities[i]
    end

    Fi = SVector{6}(tip_wrench)
    for i = n:-1:1
        Gi = SMatrix{6,6}(spatial_inertias[i])
        # The -ad(V)'*G*V term is Featherstone's v ×* (Iv): the force cross product
        # is the negative transpose of the motion cross product (twist/wrench duality).
        Fi = AdTi[i+1]' * Fi + Gi * Vdi[i+1] - ad(Vi[i+1])' * Gi * Vi[i+1]
        joint_torques[i] = Fi' * Ai[i]
    end

    return joint_torques
end

"""
    mass_matrix_crba(joint_positions, link_frames, spatial_inertias, screw_axes)

Computes the mass matrix ``M(\\theta)`` of an open chain robot using the
**Composite Rigid Body Algorithm (CRBA)**.

!!! info "Composite Rigid Body Algorithm"
    The mass matrix ``M(\\theta)`` relates joint accelerations to joint torques in the
    absence of velocity-dependent and gravitational forces:
    ``\\tau = M(\\theta) \\ddot{\\theta}``. CRBA computes ``M`` in a single backward
    pass by accumulating *composite spatial inertias* — the total inertia of each
    link plus all its descendants — from leaf to root, then assembling the matrix
    entries by projecting these composite inertias onto the joint screw axes.

!!! details "Educational note"
    The textbook (Chapter 8.2) teaches a simpler approach that calls
    [`inverse_dynamics_rnea`](@ref) ``n`` times with unit accelerations: "what torques
    are needed for unit acceleration of joint ``i``?" Each call produces one column
    of ``M(\\theta)``. This builds intuition but requires ``n`` full Newton-Euler
    passes. CRBA computes the same result with one forward pass and one backward
    pass, exploiting the symmetry of ``M``. The textbook algorithm is preserved
    See [`mass_matrix_rnea`](@ref) for this textbook algorithm.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}``.
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.

# Returns
The ``n×n`` mass matrix ``M(\\theta)``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> mass_matrix_crba([0.1, 0.1, 0.1], link_frames, spatial_inertias, screw_axes)
3×3 Matrix{Float64}:
 22.5433      -0.307147  -0.00718426
 -0.307147     1.96851    0.432157
 -0.00718426   0.432157   0.191631
```
"""
function mass_matrix_crba(
    joint_positions::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    M = zeros(n, n)
    Ai = Vector{SVector{6,Float64}}(undef, n)
    AdTi = Vector{SMatrix{6,6,Float64,36}}(undef, n + 1)
    Gc = Vector{SMatrix{6,6,Float64,36}}(undef, n)
    mass_matrix_crba!(
        M,
        joint_positions,
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Gc,
    )
end

"""
    mass_matrix_crba!(M, joint_positions, link_frames, spatial_inertias, screw_axes, Ai, AdTi, Gc)

In-place version of [`mass_matrix_crba`](@ref) using the Composite Rigid Body Algorithm (CRBA).
Writes the result into `M`. The workspace vectors `Ai` (length n), `AdTi` (length n+1),
and `Gc` (length n) must be pre-allocated.
"""
function mass_matrix_crba!(
    M::AbstractMatrix,
    joint_positions::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
    Ai::AbstractVector{SVector{6,Float64}},
    AdTi::AbstractVector{SMatrix{6,6,Float64,36}},
    Gc::AbstractVector{SMatrix{6,6,Float64,36}},
)
    n = length(joint_positions)

    # Forward pass: compute screw axes in link frames (Ai) and inter-link adjoints (AdTi)
    Mi = SMatrix{4,4,Float64}(LA.I)
    AdTi[n+1] = adjoint_representation(transform_inv(link_frames[n+1]))

    for i = 1:n
        Mi = Mi * SMatrix{4,4}(link_frames[i])
        Ai[i] =
            adjoint_representation(transform_inv(Mi)) * SVector{6}(@view screw_axes[:, i])
        AdTi[i] = adjoint_representation(
            matrix_exp6(vec_to_se3(Ai[i] * -joint_positions[i])) *
            transform_inv(SMatrix{4,4}(link_frames[i])),
        )
    end

    # Backward pass: accumulate composite spatial inertias (leaf to root)
    Gc[n] = SMatrix{6,6}(spatial_inertias[n])
    for i = (n-1):-1:1
        Gc[i] = SMatrix{6,6}(spatial_inertias[i]) + AdTi[i+1]' * Gc[i+1] * AdTi[i+1]
    end

    # Assembly: propagate F = Gc[j] * A[j] backwards through AdTi to compute
    # M[i,j] = A_i^T * F. Exploits symmetry: M[j,i] = M[i,j].
    for j = 1:n
        F = Gc[j] * Ai[j]
        M[j, j] = Ai[j]' * F
        for i = (j-1):-1:1
            F = AdTi[i+1]' * F
            M[i, j] = Ai[i]' * F
            M[j, i] = M[i, j]
        end
    end

    return M
end

"""
    mass_matrix_rnea(joint_positions, link_frames, spatial_inertias, screw_axes)

Computes the mass matrix using the textbook algorithm (Chapter 8.2): calls
[`inverse_dynamics_rnea`](@ref) ``n`` times, each with a unit acceleration vector.
This is slower than [`mass_matrix_crba`](@ref) (which uses CRBA) but directly
illustrates that column ``i`` of ``M(\\theta)`` is the torque needed for unit
acceleration of joint ``i``.
"""
function mass_matrix_rnea(
    joint_positions::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    M = zeros(n, n)
    ddq = zeros(n)
    zero_dq = zeros(n)
    zero_g = zeros(3)
    zero_F = zeros(6)
    for i = 1:n
        ddq[i] = 1
        M[:, i] = inverse_dynamics_rnea(
            joint_positions,
            zero_dq,
            ddq,
            zero_g,
            zero_F,
            link_frames,
            spatial_inertias,
            screw_axes,
        )
        ddq[i] = 0
    end
    return M
end

"""
    velocity_quadratic_forces(joint_positions, joint_velocities, link_frames, spatial_inertias, screw_axes)

Computes the Coriolis and centripetal terms ``c(\\theta, \\dot{\\theta})`` in the inverse
dynamics of an open chain robot. Calls [`inverse_dynamics_rnea`](@ref) with `gravity = 0`,
`tip_wrench = 0`, and `joint_accelerations = 0`.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `joint_velocities`: the ``n``-vector of joint rates.
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}``.
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.

# Returns
The ``n``-vector of Coriolis and centripetal terms.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> velocity_quadratic_forces([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], link_frames, spatial_inertias, screw_axes)
3-element Vector{Float64}:
  0.26453118054501235
 -0.0550515682891655
 -0.006891320068248912
```
"""
function velocity_quadratic_forces(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    tau = zeros(n)
    Ai = Vector{SVector{6,Float64}}(undef, n)
    AdTi = Vector{SMatrix{6,6,Float64,36}}(undef, n + 1)
    Vi = Vector{SVector{6,Float64}}(undef, n + 1)
    Vdi = Vector{SVector{6,Float64}}(undef, n + 1)
    velocity_quadratic_forces!(
        tau,
        joint_positions,
        joint_velocities,
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Vi,
        Vdi,
    )
end

"""
    velocity_quadratic_forces!(tau, joint_positions, joint_velocities, link_frames, spatial_inertias, screw_axes, Ai, AdTi, Vi, Vdi)

In-place version of [`velocity_quadratic_forces`](@ref). Writes the result into `tau`.
"""
function velocity_quadratic_forces!(
    tau::AbstractVector,
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
    Ai::AbstractVector{SVector{6,Float64}},
    AdTi::AbstractVector{SMatrix{6,6,Float64,36}},
    Vi::AbstractVector{SVector{6,Float64}},
    Vdi::AbstractVector{SVector{6,Float64}},
)
    # Reuse tau as temporary zero vector, then overwrite with result
    tau .= 0
    inverse_dynamics_rnea!(
        tau,
        joint_positions,
        joint_velocities,
        tau,
        SA[0.0, 0.0, 0.0],
        SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Vi,
        Vdi,
    )
end

"""
    gravity_forces(joint_positions, gravity, link_frames, spatial_inertias, screw_axes)

Computes the joint forces/torques an open chain robot requires to overcome gravity at
its configuration. Calls [`inverse_dynamics_rnea`](@ref) with `joint_velocities = 0`,
`joint_accelerations = 0`, and `tip_wrench = 0`.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `gravity`: the gravity vector ``g`` (e.g., `[0, 0, -9.8]`).
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}``.
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.

# Returns
The ``n``-vector of joint gravity torques ``g(\\theta)``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> gravity_forces([0.1, 0.1, 0.1], [0, 0, -9.8], link_frames, spatial_inertias, screw_axes)
3-element Vector{Float64}:
  28.40331261821983
 -37.64094817177068
  -5.4415891999683605
```
"""
function gravity_forces(
    joint_positions::AbstractVector,
    gravity::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    tau = zeros(n)
    Ai = Vector{SVector{6,Float64}}(undef, n)
    AdTi = Vector{SMatrix{6,6,Float64,36}}(undef, n + 1)
    Vi = Vector{SVector{6,Float64}}(undef, n + 1)
    Vdi = Vector{SVector{6,Float64}}(undef, n + 1)
    gravity_forces!(
        tau,
        joint_positions,
        gravity,
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Vi,
        Vdi,
    )
end

"""
    gravity_forces!(tau, joint_positions, gravity, link_frames, spatial_inertias, screw_axes, Ai, AdTi, Vi, Vdi)

In-place version of [`gravity_forces`](@ref). Writes the result into `tau`.
"""
function gravity_forces!(
    tau::AbstractVector,
    joint_positions::AbstractVector,
    gravity::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
    Ai::AbstractVector{SVector{6,Float64}},
    AdTi::AbstractVector{SMatrix{6,6,Float64,36}},
    Vi::AbstractVector{SVector{6,Float64}},
    Vdi::AbstractVector{SVector{6,Float64}},
)
    tau .= 0
    inverse_dynamics_rnea!(
        tau,
        joint_positions,
        tau,
        tau,
        gravity,
        SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Vi,
        Vdi,
    )
end

"""
    end_effector_forces(joint_positions, tip_wrench, link_frames, spatial_inertias, screw_axes)

Computes the joint forces/torques an open chain robot requires only to create the end-effector force `tip_wrench`.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `tip_wrench`: the spatial force applied by the end-effector expressed in frame `{n+1}`.
- `link_frames`: the list of link frames `i` relative to `i-1` at the home position.
- `spatial_inertias`: the spatial inertia matrices `Gi` of the links.
- `screw_axes`: the screw axes `Si` of the joints in a space frame, in the format of a matrix with axes as the columns.

Returns the joint forces and torques required only to create the end-effector force `tip_wrench`.
This function calls inverse_dynamics_rnea with `gravity = 0`, `joint_velocities = 0`, and `joint_accelerations = 0`.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> import LinearAlgebra as LA

julia> joint_positions = [0.1, 0.1, 0.1]
3-element Vector{Float64}:
 0.1
 0.1
 0.1

julia> tip_wrench = [1, 1, 1, 1, 1, 1]
6-element Vector{Int64}:
 1
 1
 1
 1
 1
 1

julia> M01 = [1 0 0        0;
              0 1 0        0;
              0 0 1 0.089159;
              0 0 0        1]
4×4 Matrix{Float64}:
 1.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.089159
 0.0  0.0  0.0  1.0

julia> M12 = [ 0 0 1    0.28;
               0 1 0 0.13585;
              -1 0 0       0;
               0 0 0       1]
4×4 Matrix{Float64}:
  0.0  0.0  1.0  0.28
  0.0  1.0  0.0  0.13585
 -1.0  0.0  0.0  0.0
  0.0  0.0  0.0  1.0

julia> M23 = [1 0 0       0;
              0 1 0 -0.1197;
              0 0 1   0.395;
              0 0 0       1]
4×4 Matrix{Float64}:
 1.0  0.0  0.0   0.0
 0.0  1.0  0.0  -0.1197
 0.0  0.0  1.0   0.395
 0.0  0.0  0.0   1.0

julia> M34 = [1 0 0       0;
              0 1 0       0;
              0 0 1 0.14225;
              0 0 0       1]
4×4 Matrix{Float64}:
 1.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.14225
 0.0  0.0  0.0  1.0

julia> link_frames = [M01, M12, M23, M34]
4-element Vector{Matrix{Float64}}:
 [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.089159; 0.0 0.0 0.0 1.0]
 [0.0 0.0 1.0 0.28; 0.0 1.0 0.0 0.13585; -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0]
 [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 -0.1197; 0.0 0.0 1.0 0.395; 0.0 0.0 0.0 1.0]
 [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.14225; 0.0 0.0 0.0 1.0]

julia> G1 = LA.Diagonal([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
6×6 LinearAlgebra.Diagonal{Float64, Vector{Float64}}:
 0.010267   ⋅         ⋅        ⋅    ⋅    ⋅ 
  ⋅        0.010267   ⋅        ⋅    ⋅    ⋅ 
  ⋅         ⋅        0.00666   ⋅    ⋅    ⋅ 
  ⋅         ⋅         ⋅       3.7   ⋅    ⋅ 
  ⋅         ⋅         ⋅        ⋅   3.7   ⋅ 
  ⋅         ⋅         ⋅        ⋅    ⋅   3.7

julia> G2 = LA.Diagonal([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
6×6 LinearAlgebra.Diagonal{Float64, Vector{Float64}}:
 0.22689   ⋅        ⋅          ⋅      ⋅      ⋅ 
  ⋅       0.22689   ⋅          ⋅      ⋅      ⋅ 
  ⋅        ⋅       0.0151074   ⋅      ⋅      ⋅ 
  ⋅        ⋅        ⋅         8.393   ⋅      ⋅ 
  ⋅        ⋅        ⋅          ⋅     8.393   ⋅ 
  ⋅        ⋅        ⋅          ⋅      ⋅     8.393

julia> G3 = LA.Diagonal([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])
6×6 LinearAlgebra.Diagonal{Float64, Vector{Float64}}:
 0.0494433   ⋅          ⋅         ⋅      ⋅      ⋅ 
  ⋅         0.0494433   ⋅         ⋅      ⋅      ⋅ 
  ⋅          ⋅         0.004095   ⋅      ⋅      ⋅ 
  ⋅          ⋅          ⋅        2.275   ⋅      ⋅ 
  ⋅          ⋅          ⋅         ⋅     2.275   ⋅ 
  ⋅          ⋅          ⋅         ⋅      ⋅     2.275

julia> spatial_inertias = [G1, G2, G3]
3-element Vector{LinearAlgebra.Diagonal{Float64, Vector{Float64}}}:
 Diagonal([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
 Diagonal([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
 Diagonal([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])

julia> screw_axes = [ 1  0  1      0  1      0;
                      0  1  0 -0.089  0      0;
                      0  1  0 -0.089  0  0.425]'
6×3 adjoint(::Matrix{Float64}) with eltype Float64:
 1.0   0.0     0.0
 0.0   1.0     1.0
 1.0   0.0     0.0
 0.0  -0.089  -0.089
 1.0   0.0     0.0
 0.0   0.0     0.425

julia> end_effector_forces(joint_positions, tip_wrench, link_frames, spatial_inertias, screw_axes)
3-element Vector{Float64}:
 1.4095460782639782
 1.8577149723180628
 1.392409
```
"""
function end_effector_forces(
    joint_positions::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    tau = zeros(n)
    Ai = Vector{SVector{6,Float64}}(undef, n)
    AdTi = Vector{SMatrix{6,6,Float64,36}}(undef, n + 1)
    Vi = Vector{SVector{6,Float64}}(undef, n + 1)
    Vdi = Vector{SVector{6,Float64}}(undef, n + 1)
    end_effector_forces!(
        tau,
        joint_positions,
        tip_wrench,
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Vi,
        Vdi,
    )
end

"""
    end_effector_forces!(tau, joint_positions, tip_wrench, link_frames, spatial_inertias, screw_axes, Ai, AdTi, Vi, Vdi)

In-place version of [`end_effector_forces`](@ref). Writes the result into `tau`.
"""
function end_effector_forces!(
    tau::AbstractVector,
    joint_positions::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
    Ai::AbstractVector{SVector{6,Float64}},
    AdTi::AbstractVector{SMatrix{6,6,Float64,36}},
    Vi::AbstractVector{SVector{6,Float64}},
    Vdi::AbstractVector{SVector{6,Float64}},
)
    tau .= 0
    inverse_dynamics_rnea!(
        tau,
        joint_positions,
        tau,
        tau,
        SA[0.0, 0.0, 0.0],
        tip_wrench,
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Vi,
        Vdi,
    )
end

"""
    forward_dynamics_crba(joint_positions, joint_velocities, joint_torques, gravity, tip_wrench, link_frames, spatial_inertias, screw_axes)

Computes forward dynamics in the space frame for an open chain robot.

!!! info "How does forward dynamics work?"
    Forward dynamics answers: "given the applied torques, how does the robot accelerate?"
    It solves the manipulator equation for joint accelerations:

    ```math
    \\ddot{\\theta} = M(\\theta)^{-1} \\left[ \\tau - c(\\theta,\\dot{\\theta}) - g(\\theta) - J^T(\\theta) \\mathcal{F}_{\\text{tip}} \\right]
    ```

    where ``M`` is the [mass matrix](@ref mass_matrix_crba) (computed via CRBA), and the
    right-hand side is computed via a single [`inverse_dynamics_rnea`](@ref) call with zero
    accelerations that captures velocity-dependent, gravitational, and tip wrench forces.

!!! details "Educational note"
    The textbook computes forward dynamics by explicitly forming ``M^{-1}`` and calling
    [`inverse_dynamics_rnea`](@ref) separately for Coriolis, gravity, and tip wrench terms.
    This implementation makes two improvements:

    1. **Single RNEA call**: instead of computing ``c``, ``g``, and ``J^T F`` separately,
       a single call to RNEA with zero accelerations returns ``\\tau_{\\text{bias}} = c + g + J^T F``.
    2. **Backslash solve** (`M \\ b`): instead of explicitly forming ``M^{-1}`` and
       multiplying, Julia's `\\` operator solves the linear system ``M x = b`` via LU
       factorization. This is both faster (avoids forming the inverse) and more
       numerically stable (fewer floating-point operations, better conditioning).

    For the fastest forward dynamics, use [`forward_dynamics_aba`](@ref) which avoids
    forming ``M`` entirely. See [`forward_dynamics_rnea`](@ref) for the textbook algorithm.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `joint_velocities`: the ``n``-vector of joint rates.
- `joint_torques`: the ``n``-vector of joint forces/torques.
- `gravity`: the gravity vector ``g`` (e.g., `[0, 0, -9.8]`).
- `tip_wrench`: the wrench ``\\mathcal{F}_{\\text{tip}}`` applied by the end-effector expressed in frame ``\\{n+1\\}``.
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}``.
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.

# Returns
The ``n``-vector of joint accelerations ``\\ddot{\\theta}``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> forward_dynamics_crba([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], [0.5, 0.6, 0.7], [0, 0, -9.8], [1, 1, 1, 1, 1, 1], link_frames, spatial_inertias, screw_axes)
3-element Vector{Float64}:
  -0.9739290670855625
  25.584667840340547
 -32.91499212478147
```
"""
function forward_dynamics_crba(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_torques::AbstractVector,
    gravity::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    ddq = zeros(n)
    M = zeros(n, n)
    tau_bias = zeros(n)
    Ai = Vector{SVector{6,Float64}}(undef, n)
    AdTi = Vector{SMatrix{6,6,Float64,36}}(undef, n + 1)
    Gc = Vector{SMatrix{6,6,Float64,36}}(undef, n)
    Vi = Vector{SVector{6,Float64}}(undef, n + 1)
    Vdi = Vector{SVector{6,Float64}}(undef, n + 1)
    forward_dynamics_crba!(
        ddq,
        joint_positions,
        joint_velocities,
        joint_torques,
        gravity,
        tip_wrench,
        link_frames,
        spatial_inertias,
        screw_axes,
        M,
        tau_bias,
        Ai,
        AdTi,
        Gc,
        Vi,
        Vdi,
    )
end

"""
    forward_dynamics_crba!(ddq, joint_positions, joint_velocities, joint_torques, gravity, tip_wrench, link_frames, spatial_inertias, screw_axes, M, tau_bias, Ai, AdTi, Gc, Vi, Vdi)

In-place version of [`forward_dynamics_crba`](@ref). Writes the result into `ddq`.
Workspace: `M` (n×n), `tau_bias` (n), `Ai` (n), `AdTi` (n+1), `Gc` (n), `Vi` (n+1), `Vdi` (n+1).
"""
function forward_dynamics_crba!(
    ddq::AbstractVector,
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_torques::AbstractVector,
    gravity::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
    M::AbstractMatrix,
    tau_bias::AbstractVector,
    Ai::AbstractVector{SVector{6,Float64}},
    AdTi::AbstractVector{SMatrix{6,6,Float64,36}},
    Gc::AbstractVector{SMatrix{6,6,Float64,36}},
    Vi::AbstractVector{SVector{6,Float64}},
    Vdi::AbstractVector{SVector{6,Float64}},
)
    mass_matrix_crba!(
        M,
        joint_positions,
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Gc,
    )
    # Reuse ddq as temporary zero-acceleration vector, then overwrite with result
    ddq .= 0
    inverse_dynamics_rnea!(
        tau_bias,
        joint_positions,
        joint_velocities,
        ddq,
        gravity,
        tip_wrench,
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Vi,
        Vdi,
    )
    ddq .= M \ (joint_torques - tau_bias)
    return ddq
end

"""
    forward_dynamics_aba(joint_positions, joint_velocities, joint_torques, gravity, tip_wrench, link_frames, spatial_inertias, screw_axes)

Computes forward dynamics using the **Articulated Body Algorithm (ABA)**.

!!! info "Articulated Body Algorithm"
    ABA solves for joint accelerations in O(n) with three passes over the kinematic
    chain, without forming or inverting the mass matrix:

    1. **Outward pass**: compute link velocities and velocity-dependent bias forces.
    2. **Inward pass**: accumulate *articulated body inertias* — the effective inertia
       of each subtree accounting for joint freedom — from leaf to root.
    3. **Outward pass**: solve for joint accelerations from root to leaf.

    This is the same algorithm Pinocchio uses for forward dynamics.

!!! details "Educational note"
    The textbook (Chapter 8.3) computes forward dynamics by explicitly forming
    ``M(\\theta)`` and solving ``M \\ddot{\\theta} = \\tau - c - g``. See
    [`forward_dynamics_crba`](@ref) for that approach. ABA avoids forming ``M``
    entirely, which is asymptotically faster for large ``n``.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `joint_velocities`: the ``n``-vector of joint rates.
- `joint_torques`: the ``n``-vector of joint forces/torques.
- `gravity`: the gravity vector ``g`` (e.g., `[0, 0, -9.8]`).
- `tip_wrench`: the wrench ``\\mathcal{F}_{\\text{tip}}`` applied by the end-effector expressed in frame ``\\{n+1\\}``.
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}``.
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.

# Returns
The ``n``-vector of joint accelerations ``\\ddot{\\theta}``.
"""
function forward_dynamics_aba(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_torques::AbstractVector,
    gravity::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    ddq = zeros(n)
    Ai = Vector{SVector{6,Float64}}(undef, n)
    AdTi = Vector{SMatrix{6,6,Float64,36}}(undef, n + 1)
    Vi = Vector{SVector{6,Float64}}(undef, n + 1)
    ci = Vector{SVector{6,Float64}}(undef, n)
    IA = Vector{SMatrix{6,6,Float64,36}}(undef, n)
    pA = Vector{SVector{6,Float64}}(undef, n)
    U = Vector{SVector{6,Float64}}(undef, n)
    D = Vector{Float64}(undef, n)
    u = Vector{Float64}(undef, n)
    forward_dynamics_aba!(
        ddq,
        joint_positions,
        joint_velocities,
        joint_torques,
        gravity,
        tip_wrench,
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Vi,
        ci,
        IA,
        pA,
        U,
        D,
        u,
    )
end

"""
    forward_dynamics_aba!(ddq, joint_positions, joint_velocities, joint_torques, gravity, tip_wrench, link_frames, spatial_inertias, screw_axes, Ai, AdTi, Vi, ci, IA, pA, U, D, u)

In-place version of [`forward_dynamics_aba`](@ref). Writes the result into `ddq`.
The workspace vectors must be pre-allocated: `Ai` (n), `AdTi` (n+1), `Vi` (n+1),
`ci` (n), `IA` (n), `pA` (n), `U` (n), `D` (n), `u` (n).
"""
function forward_dynamics_aba!(
    ddq::AbstractVector,
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_torques::AbstractVector,
    gravity::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
    Ai::AbstractVector{SVector{6,Float64}},
    AdTi::AbstractVector{SMatrix{6,6,Float64,36}},
    Vi::AbstractVector{SVector{6,Float64}},
    ci::AbstractVector{SVector{6,Float64}},
    IA::AbstractVector{SMatrix{6,6,Float64,36}},
    pA::AbstractVector{SVector{6,Float64}},
    U::AbstractVector{SVector{6,Float64}},
    D::AbstractVector{Float64},
    u::AbstractVector{Float64},
)
    n = length(joint_positions)

    # Pass 1 (outward): velocities, bias accelerations, bias forces
    Mi = SMatrix{4,4,Float64}(LA.I)
    Vi[1] = @SVector zeros(6)
    AdTi[n+1] = adjoint_representation(transform_inv(link_frames[n+1]))

    for i = 1:n
        Mi = Mi * SMatrix{4,4}(link_frames[i])
        Ai[i] =
            adjoint_representation(transform_inv(Mi)) * SVector{6}(@view screw_axes[:, i])
        AdTi[i] = adjoint_representation(
            matrix_exp6(vec_to_se3(Ai[i] * -joint_positions[i])) *
            transform_inv(SMatrix{4,4}(link_frames[i])),
        )
        Vi[i+1] = AdTi[i] * Vi[i] + Ai[i] * joint_velocities[i]
        ci[i] = ad(Vi[i+1]) * Ai[i] * joint_velocities[i]
        Gi = SMatrix{6,6}(spatial_inertias[i])
        IA[i] = Gi
        # Featherstone's velocity-product force is p_A = v ×* (I v), where ×* is the
        # spatial force cross product. Twists and wrenches live in dual spaces, so the
        # force cross product is the negative transpose of the motion cross product:
        # v ×* = -ad(V)'. Hence pA = v ×* (G V) = -ad(V)' * G * V.
        pA[i] = -ad(Vi[i+1])' * Gi * Vi[i+1]
    end

    pA[n] = pA[n] + AdTi[n+1]' * SVector{6}(tip_wrench)

    # Pass 2 (inward): articulated body inertias and bias forces
    for i = n:-1:1
        U[i] = IA[i] * Ai[i]
        D[i] = Ai[i]' * U[i]
        u[i] = joint_torques[i] - Ai[i]' * pA[i]
        if i > 1
            Ia = IA[i] - U[i] * U[i]' / D[i]
            pa = pA[i] + Ia * ci[i] + U[i] * u[i] / D[i]
            IA[i-1] = IA[i-1] + AdTi[i]' * Ia * AdTi[i]
            pA[i-1] = pA[i-1] + AdTi[i]' * pa
        end
    end

    # Pass 3 (outward): joint accelerations
    a0 = SA[0.0, 0.0, 0.0, -gravity[1], -gravity[2], -gravity[3]]
    a = AdTi[1] * a0

    for i = 1:n
        a_plus_c = a + ci[i]
        ddq[i] = (u[i] - U[i]' * a_plus_c) / D[i]
        if i < n
            a = AdTi[i+1] * (a_plus_c + Ai[i] * ddq[i])
        end
    end

    return ddq
end

"""
    forward_dynamics_rnea(joint_positions, joint_velocities, joint_torques, gravity, tip_wrench, link_frames, spatial_inertias, screw_axes)

Computes forward dynamics using the textbook algorithm (Chapter 8.3): explicitly
forms ``M^{-1}`` and calls [`inverse_dynamics_rnea`](@ref) separately for Coriolis,
gravity, and tip wrench terms. This is slower than [`forward_dynamics_crba`](@ref)
(which uses CRBA + single RNEA + backslash) but directly mirrors the textbook
equation ``\\ddot{\\theta} = M^{-1}[\\tau - c - g - J^T F_{\\text{tip}}]``.
"""
function forward_dynamics_rnea(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_torques::AbstractVector,
    gravity::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    n = length(joint_positions)
    ddq = zeros(n)
    M = zeros(n, n)
    tau_c = zeros(n)
    tau_g = zeros(n)
    tau_f = zeros(n)
    Ai = Vector{SVector{6,Float64}}(undef, n)
    AdTi = Vector{SMatrix{6,6,Float64,36}}(undef, n + 1)
    Gc = Vector{SMatrix{6,6,Float64,36}}(undef, n)
    Vi = Vector{SVector{6,Float64}}(undef, n + 1)
    Vdi = Vector{SVector{6,Float64}}(undef, n + 1)
    forward_dynamics_rnea!(
        ddq,
        joint_positions,
        joint_velocities,
        joint_torques,
        gravity,
        tip_wrench,
        link_frames,
        spatial_inertias,
        screw_axes,
        M,
        tau_c,
        tau_g,
        tau_f,
        Ai,
        AdTi,
        Gc,
        Vi,
        Vdi,
    )
end

"""
    forward_dynamics_rnea!(ddq, joint_positions, joint_velocities, joint_torques, gravity, tip_wrench, link_frames, spatial_inertias, screw_axes, M, tau_c, tau_g, tau_f, Ai, AdTi, Gc, Vi, Vdi)

In-place version of [`forward_dynamics_rnea`](@ref). Writes the result into `ddq`.
Workspace: `M` (n×n), `tau_c` (n), `tau_g` (n), `tau_f` (n), `Ai` (n), `AdTi` (n+1),
`Gc` (n), `Vi` (n+1), `Vdi` (n+1).
"""
function forward_dynamics_rnea!(
    ddq::AbstractVector,
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_torques::AbstractVector,
    gravity::AbstractVector,
    tip_wrench::AbstractVector,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
    M::AbstractMatrix,
    tau_c::AbstractVector,
    tau_g::AbstractVector,
    tau_f::AbstractVector,
    Ai::AbstractVector{SVector{6,Float64}},
    AdTi::AbstractVector{SMatrix{6,6,Float64,36}},
    Gc::AbstractVector{SMatrix{6,6,Float64,36}},
    Vi::AbstractVector{SVector{6,Float64}},
    Vdi::AbstractVector{SVector{6,Float64}},
)
    zero_n = ddq
    zero_3 = SA[0.0, 0.0, 0.0]
    zero_6 = SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # Mass matrix via CRBA
    mass_matrix_crba!(
        M,
        joint_positions,
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Gc,
    )

    # Coriolis/centrifugal: RNEA with gravity=0, tip_wrench=0, ddq=0
    zero_n .= 0
    inverse_dynamics_rnea!(
        tau_c,
        joint_positions,
        joint_velocities,
        zero_n,
        zero_3,
        zero_6,
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Vi,
        Vdi,
    )

    # Gravity: RNEA with dq=0, ddq=0, tip_wrench=0
    zero_n .= 0
    inverse_dynamics_rnea!(
        tau_g,
        joint_positions,
        zero_n,
        zero_n,
        gravity,
        zero_6,
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Vi,
        Vdi,
    )

    # End-effector forces: RNEA with dq=0, ddq=0, gravity=0
    zero_n .= 0
    inverse_dynamics_rnea!(
        tau_f,
        joint_positions,
        zero_n,
        zero_n,
        zero_3,
        tip_wrench,
        link_frames,
        spatial_inertias,
        screw_axes,
        Ai,
        AdTi,
        Vi,
        Vdi,
    )

    # ddq = M⁻¹ * (τ - c - g - f)
    ddq .= joint_torques .- tau_c .- tau_g .- tau_f
    ddq .= LA.inv(M) * ddq
    return ddq
end

"""
    euler_step(joint_positions, joint_velocities, joint_accelerations, timestep)

Compute the joint angles and velocities at the next timestep using first order Euler integration.

# Arguments
- `joint_positions`: the ``n``-vector of joint variables.
- `joint_velocities`: the ``n``-vector of joint rates.
- `joint_accelerations`: the ``n``-vector of joint accelerations.
- `timestep`: the timestep delta t.

# Return
- `joint_positionsNext`: the vector of joint variables after `timestep` from first order Euler integration.
- `joint_velocitiesNext`: the vector of joint rates after `timestep` from first order Euler integration.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> euler_step([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], [2, 1.5, 1], 0.1)
([0.11000000000000001, 0.12000000000000001, 0.13], [0.30000000000000004, 0.35000000000000003, 0.4])
```
"""
function euler_step(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_accelerations::AbstractVector,
    timestep::Number,
)
    joint_positions + timestep * joint_velocities,
    joint_velocities + timestep * joint_accelerations
end

"""
    inverse_dynamics_trajectory(joint_position_traj, joint_velocity_traj, joint_acceleration_traj, gravity, tip_wrench_traj, link_frames, spatial_inertias, screw_axes)

Calculates the joint forces/torques required to move the serial chain along the given
trajectory using [`inverse_dynamics_rnea`](@ref) at each timestep.

# Arguments
- `joint_position_traj`: an ``N×n`` matrix of joint variables, where row ``i`` is the joint vector at timestep ``i``.
- `joint_velocity_traj`: an ``N×n`` matrix of joint velocities.
- `joint_acceleration_traj`: an ``N×n`` matrix of joint accelerations.
- `gravity`: the gravity vector ``g`` (e.g., `[0, 0, -9.8]`).
- `tip_wrench_traj`: an ``N×6`` matrix where row ``i`` is the spatial wrench applied by the end-effector at timestep ``i``.
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}``.
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.

# Returns
An ``N×n`` matrix of joint forces/torques, where row ``i`` is the joint force-torque at timestep ``i``.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> traj_pos = [0.0 0.0 0.0; 0.5 0.5 0.5; 1.0 1.0 1.0];

julia> traj_vel = [0.0 0.0 0.0; 0.5 0.5 0.5; 1.0 1.0 1.0];

julia> traj_acc = [0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5];

julia> inverse_dynamics_trajectory(traj_pos, traj_vel, traj_acc, [0, 0, -9.8], ones(3, 6), link_frames, spatial_inertias, screw_axes)
3×3 adjoint(::Matrix{Float64}) with eltype Float64:
  23.7113  -35.2321  -3.87345
 105.336   -27.6441  -1.41169
 129.143   -17.3084   2.30991
```
"""
function inverse_dynamics_trajectory(
    joint_position_traj::AbstractMatrix,
    joint_velocity_traj::AbstractMatrix,
    joint_acceleration_traj::AbstractMatrix,
    gravity::AbstractVector,
    tip_wrench_traj::AbstractMatrix,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
)
    joint_position_traj = joint_position_traj'
    joint_velocity_traj = joint_velocity_traj'
    joint_acceleration_traj = joint_acceleration_traj'
    tip_wrench_traj = tip_wrench_traj'
    joint_torque_traj = copy(joint_position_traj)

    for i in axes(joint_position_traj, 2)
        joint_torque_traj[:, i] = inverse_dynamics_rnea(
            joint_position_traj[:, i],
            joint_velocity_traj[:, i],
            joint_acceleration_traj[:, i],
            gravity,
            tip_wrench_traj[:, i],
            link_frames,
            spatial_inertias,
            screw_axes,
        )
    end

    joint_torque_traj'
end

"""
    forward_dynamics_trajectory(joint_positions, joint_velocities, joint_torque_traj, gravity, tip_wrench_traj, link_frames, spatial_inertias, screw_axes, timestep, integration_resolution)

Simulates the motion of a serial chain given an open-loop history of joint
forces/torques. Uses [`forward_dynamics_crba`](@ref) with [`euler_step`](@ref) integration at
each timestep.

# Arguments
- `joint_positions`: the ``n``-vector of initial joint variables.
- `joint_velocities`: the ``n``-vector of initial joint velocities.
- `joint_torque_traj`: an ``N×n`` matrix where row ``i`` is the joint force/torque vector at timestep ``i``.
- `gravity`: the gravity vector ``g`` (e.g., `[0, 0, -9.8]`).
- `tip_wrench_traj`: an ``N×6`` matrix where row ``i`` is the spatial wrench applied by the end-effector at timestep ``i``.
- `link_frames`: a vector of ``n+1`` SE(3) matrices, where `link_frames[i]` is ``M_{i-1,i}`` and `link_frames[n+1]` is ``M_{n,n+1}``.
- `spatial_inertias`: a vector of ``n`` symmetric 6×6 spatial inertia matrices ``G_i`` of the links.
- `screw_axes`: the screw axes ``S_i`` of the joints in a space frame, as a 6×``n`` matrix with axes as columns.
- `timestep`: the timestep ``\\Delta t`` between consecutive trajectory points.
- `integration_resolution`: the number of Euler integration steps during each timestep ``\\Delta t``. Must be a positive integer ``\\geq 1``. Larger values result in slower simulations but less accumulation of integration error.

# Returns
- `joint_position_traj`: the ``N×n`` matrix of joint variables resulting from the specified joint forces/torques.
- `joint_velocity_traj`: the ``N×n`` matrix of joint velocities resulting from the specified joint forces/torques.

# Examples
```jldoctest; setup = :(using ModernRoboticsBook)
julia> joint_torque_traj = [3.63 -6.58 -5.57; 3.74 -5.55 -5.5; 4.31 -0.68 -5.19];

julia> thetamat, dthetamat = forward_dynamics_trajectory([0.1, 0.1, 0.1], [0.1, 0.2, 0.3], joint_torque_traj, [0, 0, -9.8], ones(3, 6), link_frames, spatial_inertias, screw_axes, 0.1, 8);

julia> thetamat
3×3 adjoint(::Matrix{Float64}) with eltype Float64:
 0.1       0.1        0.1
 0.106431  0.2626    -0.226649
 0.10198   0.715813  -1.22522
```
"""
function forward_dynamics_trajectory(
    joint_positions::AbstractVector,
    joint_velocities::AbstractVector,
    joint_torque_traj::AbstractMatrix,
    gravity::AbstractVector,
    tip_wrench_traj::AbstractMatrix,
    link_frames::AbstractVector,
    spatial_inertias::AbstractVector,
    screw_axes::AbstractMatrix,
    timestep::Number,
    integration_resolution::Number,
)
    joint_torque_traj = joint_torque_traj'
    tip_wrench_traj = tip_wrench_traj'
    joint_position_traj = copy(joint_torque_traj)
    joint_position_traj[:, 1] = joint_positions
    joint_velocity_traj = copy(joint_torque_traj)
    joint_velocity_traj[:, 1] = joint_velocities

    for i = first(axes(joint_torque_traj, 2)):(last(axes(joint_torque_traj, 2))-1)
        for j = 1:integration_resolution
            joint_accelerations = forward_dynamics_crba(
                joint_positions,
                joint_velocities,
                joint_torque_traj[:, i],
                gravity,
                tip_wrench_traj[:, i],
                link_frames,
                spatial_inertias,
                screw_axes,
            )
            joint_positions, joint_velocities = euler_step(
                joint_positions,
                joint_velocities,
                joint_accelerations,
                1.0 * timestep / integration_resolution,
            )
        end

        joint_position_traj[:, i+1] = joint_positions
        joint_velocity_traj[:, i+1] = joint_velocities
    end

    joint_position_traj = joint_position_traj'
    joint_velocity_traj = joint_velocity_traj'
    return joint_position_traj, joint_velocity_traj
end

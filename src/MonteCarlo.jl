mutable struct MC <: AbstractMC
    Ham::Hamiltonian
    kappa_up::Vector{Int}
    kappa_down::Vector{Int}
    W_up::AbstractMatrix
    W_down::AbstractMatrix
    # Caches for zero-allocation updates
    W_up_col_cache::AbstractVector
    W_up_row_cache::AbstractVector
    W_down_col_cache::AbstractVector
    W_down_row_cache::AbstractVector
end

function reevaluateW!(mc::MC)
    # Calculate inverse matrices
    tilde_U_up = tilde_U(mc.Ham.U_up, mc.kappa_up)
    tilde_U_down = tilde_U(mc.Ham.U_down, mc.kappa_down)
    U_upinvs = tilde_U_up \ I
    U_downinvs = tilde_U_down \ I
    # Calculate W matrices using matrix multiplication
    mc.W_up = mc.Ham.U_up * U_upinvs
    mc.W_down = mc.Ham.U_down * U_downinvs

    return nothing
end


"""
    tilde_U(U::AbstractMatrix, kappa::Vector{Int})
------------
Creates a tilde matrix by rearranging rows of U according to kappa indices.

Parameters:
- `U`: Source matrix of size (n × m)
- `kappa`: Vector of indices where each non-zero value l indicates that row Rl of U
          should be placed at row l of the output

Returns:
- A matrix of size (m × m) with same element type as U
"""
function tilde_U(U::AbstractMatrix, kappa::Vector{Int})
    m = size(U, 2)
    n = size(U, 1)
    # check if kappa is valid
    length(kappa) == n || throw(
        DimensionMismatch(
            "Length of kappa ($(length(kappa))) must match number of rows in U ($n)",
        ),
    )
    length(filter(x -> x != 0, kappa)) == m ||
        throw(ArgumentError("kappa ($kappa) is not valid"))

    # Create output matrix with same element type as U and requested size
    tilde_U = zeros(eltype(U), m, m)

    @inbounds for (Rl, l) in enumerate(kappa)
        if l != 0
            (1 ≤ l ≤ m) || throw(BoundsError(tilde_U, (l, :)))
            tilde_U[l, :] = U[Rl, :]
        end
    end

    return tilde_U
end

"""
    is_occupied(kappa::Vector{Int}, l::Int) -> Bool

Check if site `l` is occupied in the kappa configuration vector.

Throws:
    BoundsError: if l is outside the valid range of kappa
"""
@inline function is_occupied(kappa::Vector{Int}, l::Int)
    @boundscheck 1 ≤ l ≤ length(kappa) || throw(BoundsError(kappa, l))
    @inbounds return kappa[l] != 0
end

"""
    MC(params::AbstractDict)
------------
Create a Monte Carlo object from a dictionary of parameters.
This is the user-facing, high-level constructor.
"""
function MC(params::AbstractDict)
    n1 = params[:n1]
    n2 = params[:n2]
    PBC = params[:PBC]
    antiPBC = get(params, :antiPBC, (false, false))
    lat_type = get(params, :lattice, DoubleKagome)
    B = get(params, :B, 0.0)
    lat = lat_type(1.0, n1, n2, PBC; antiPBC = antiPBC)
    N_up = params[:N_up]
    N_down = params[:N_down]
    link_in = get(params, :link_in, pi_link_in)
    link_inter = get(params, :link_inter, pi_link_inter)
    Ham = Hamiltonian(N_up, N_down, lat; link_in = link_in, link_inter = link_inter, B = B)
    ns = n1 * n2 * 3
    kappa_up = zeros(Int, ns)
    kappa_down = zeros(Int, ns)
    W_up = zeros(eltype(Ham.U_up), ns, N_up)
    W_down = zeros(eltype(Ham.U_down), ns, N_down)

    return MC(Ham, kappa_up, kappa_down, W_up, W_down)
end

"""
    MC(Ham, kappa_up, kappa_down, W_up, W_down)
------------
Create a Monte Carlo object from its core components.
This constructor automatically creates the necessary cache arrays.
It's useful for internal logic and testing.
"""
function MC(
    Ham::Hamiltonian,
    kappa_up::Vector{Int},
    kappa_down::Vector{Int},
    W_up::AbstractMatrix,
    W_down::AbstractMatrix,
)
    ns, N_up = size(W_up)
    _, N_down = size(W_down)
    @assert ns != 0
    @assert N_up != 0
    @assert N_down != 0


    return MC(
        Ham,
        kappa_up,
        kappa_down,
        W_up,
        W_down,
        zeros(eltype(W_up), ns),      # W_up_col_cache
        zeros(eltype(W_up), N_up),     # W_up_row_cache
        zeros(eltype(W_down), ns),      # W_down_col_cache
        zeros(eltype(W_down), N_down),   # W_down_row_cache
    )
end

"""
    update_W_matrices(mc::MC; K_up::Int, K_down::Int, l_up::Int, l_down::Int)
------------
Update the W matrices
"""
function update_W_matrices!(mc::MC; K_up::Int, K_down::Int, l_up::Int, l_down::Int)
    update_W!(mc.W_up, l_up, K_up, mc.W_up_col_cache, mc.W_up_row_cache)
    update_W!(mc.W_down, l_down, K_down, mc.W_down_col_cache, mc.W_down_row_cache)
end

"""
    update_W!(W::AbstractMatrix; l::Int, K::Int)
------------
Update the W matrix
``W'_{I,j} = W_{I,j} - W_{I,l} / W_{K,l} * (W_{K,j} - \\delta_{l,j})``
"""
function update_W!(
    W::AbstractMatrix,
    l::Int,
    K::Int,
    col_cache::AbstractVector,
    row_cache::AbstractVector,
)
    copyto!(row_cache, view(W, K, :))
    row_cache[l] -= 1.0
    copyto!(col_cache, view(W, :, l))
    alpha = -one(eltype(W)) / W[K, l]
    LinearAlgebra.BLAS.geru!(alpha, col_cache, row_cache, W)
    return nothing
end


"""
    init_conf_qr!(mc::MC, ns::Int, N_up::Int)

Ref: Quantum Monte Carlo Approaches for Correlated Systems (Becca and Sorella, 2017) P130

As system size increases, ``<Φ|x>`` becomes exponentially small comparing to system size, inversely proportional to the dimension of the Hilbert space.

The smallness of init ``<Φ|x>`` will lead to numerical instability in the first computation of the W matrices, but not affecting the Markov chain sampling, since we always calculate ``\\frac{<Φ|x'>}{<Φ|x>}`` there, which shall not be a numerically instable number.

To avoid this issue we initializes the particle configurations `kappa_up` and `kappa_down` in the `mc` object
using a deterministic method based on QR decomposition with column pivoting.

This approach is designed to select a set of sites for the up and down spin
electrons that corresponds to a set of linearly independent rows from the `U_up` and
`U_down` matrices, respectively. This maximizes the chances of producing a
non-singular `tilde_U` matrix, avoiding the need for random trial-and-error.

The procedure is as follows:
1.  It performs QR decomposition on `mc.Ham.U_up'` to select the `N_up` most
    linearly independent rows (sites) for the up-spin electrons.
2.  It then takes the remaining, unoccupied sites and performs a second QR
    decomposition on the corresponding subset of `mc.Ham.U_down'` to select the
    `N_down` sites for the down-spin electrons.

Parameters:
- `mc`: The Monte Carlo object, which will be modified in-place.
- `ns`: The total number of sites in the lattice.
- `N_up`: The number of up-spin electrons.
"""
function init_conf_qr!(mc::MC, ns::Int, N_up::Int)
    # Use QR decomposition with pivoting to select linearly independent rows from U
    # This should give a non-singular tilde_U matrix if possible.

    # For kappa_up
    Q_up, R_up, p_up = qr(mc.Ham.U_up', ColumnNorm())
    sites_up = p_up[1:N_up]
    mc.kappa_up = zeros(Int, ns)
    for (i, site) in enumerate(sites_up)
        mc.kappa_up[site] = i
    end

    # For kappa_down
    N_down = ns - N_up
    if N_down > 0
        available_sites = setdiff(1:ns, sites_up)
        U_down_subset = mc.Ham.U_down[available_sites, :]

        Q_down, R_down, p_down_subset = qr(U_down_subset', ColumnNorm())
        sites_down_indices = p_down_subset[1:N_down]
        sites_down = available_sites[sites_down_indices]

        mc.kappa_down = zeros(Int, ns)
        for (i, site) in enumerate(sites_down)
            mc.kappa_down[site] = i
        end
    else
        mc.kappa_down = zeros(Int, ns)
    end

    return nothing
end

"""
    find_initial_configuration!(mc::MC, rng::AbstractRNG, ns::Int, N_up::Int, n1::Int)

Finds a non-singular initial configuration for the MC simulation.
This version uses a deterministic QR-based method and should not require retries.
"""
function find_initial_configuration!(mc::MC, ns::Int, N_up::Int)
    init_conf_qr!(mc, ns, N_up)

    tilde_U_up = tilde_U(mc.Ham.U_up, mc.kappa_up)
    N_down = ns - N_up
    if N_down > 0
        tilde_U_down = tilde_U(mc.Ham.U_down, mc.kappa_down)
    end

    try
        U_upinvs = tilde_U_up \ I
        mc.W_up = mc.Ham.U_up * U_upinvs

        if N_down > 0
            tilde_U_down = tilde_U(mc.Ham.U_down, mc.kappa_down)
            U_downinvs = tilde_U_down \ I
            mc.W_down = mc.Ham.U_down * U_downinvs
        end
        return # Success
    catch e
        if e isa SingularException || e isa LinearAlgebra.LAPACKException
            error(
                "QR-based configuration is singular. The Hamiltonian may be rank-deficient.",
            )
        else
            rethrow(e)
        end
    end
end


"""
    Carlo.init!(mc::MC, ctx::MCContext, params::AbstractDict)
------------
Initialize the Monte Carlo object
`params`
* `n1` : `Int` number of cells in x direction
* `n2` : `Int` number of cells in y direction
* `PBC` : `Tuple{Bool,2}` boundary condition, e.g. (false, false)
"""
@inline function Carlo.init!(mc::MC, ctx::MCContext, params::AbstractDict)
    n1 = params[:n1]
    n2 = params[:n2]
    ns = n1 * n2 * 3
    N_up = params[:N_up]
    find_initial_configuration!(mc, ns, N_up)
end

function Z(nn::AbstractArray, kappa_up::AbstractVector, kappa_down::AbstractVector)
    # iterate over all possible moves
    # if two sites connected by a bond, check if they are occupied by different spins
    count = 0
    @inbounds for bond in nn
        site1 = bond[1]
        site2 = bond[2]
        if kappa_up[site1] != 0 && kappa_down[site2] != 0
            count += 1
        elseif kappa_up[site2] != 0 && kappa_down[site1] != 0
            count += 1
        end
    end
    return count
end

function update_configurations!(mc, flag::Int, i::Int, site::Int, l_up::Int, l_down::Int)
    if flag == 1
        # Update W matrices
        update_W_matrices!(mc; K_up = site, K_down = i, l_up = l_up, l_down = l_down)
        # Update kappa configurations
        mc.kappa_up[i], mc.kappa_up[site] = 0, l_up
        mc.kappa_down[i], mc.kappa_down[site] = l_down, 0
    else
        # Update W matrices
        update_W_matrices!(mc; K_up = i, K_down = site, l_up = l_up, l_down = l_down)
        # Update kappa configurations
        mc.kappa_up[i], mc.kappa_up[site] = l_up, 0
        mc.kappa_down[i], mc.kappa_down[site] = 0, l_down
    end
end

"""
    Carlo.sweep!(mc::MC, ctx::MCContext) -> Nothing

Perform one Monte Carlo sweep for the Mott state simulation. This implements
a two-spin swap update where electrons can exchange positions between different
spin states at occupied sites.

Note:
The W matrices are re-evaluated periodically (every n_occupied sweeps) to
maintain numerical stability. If the re-evaluation fails due to singular
matrices, a warning is issued and the simulation continues.
"""
@inline function Carlo.sweep!(mc::MC, ctx::MCContext)
    # MCMC scheme for Mott state
    # Electrons swap between different spins at occupied sites
    nn = mc.Ham.nn
    ns = length(mc.kappa_up)
    # calculate the number of neighbor states Z_{\mu}
    Zmu = Z(nn, mc.kappa_up, mc.kappa_down)
    Zmax = length(nn)
    r = rand(ctx.rng)
    if r > Zmu / Zmax
        measure!(ctx, :acc, 0.0)
        return nothing
    end
    # Randomly select neighboring sites
    bond = sample(ctx.rng, nn)
    i = bond[1]
    site = bond[2]
    ratio = 0

    # Check possible updates: can swap if both sites are occupied by different spins
    update_pool = [
        is_occupied(mc.kappa_up, i) && is_occupied(mc.kappa_down, site),
        is_occupied(mc.kappa_up, site) && is_occupied(mc.kappa_down, i),
    ]

    maybe_update = findall(identity, update_pool)
    if isempty(maybe_update)
        measure!(ctx, :acc, 0.0)
        return nothing
    end

    flag = sample(ctx.rng, maybe_update)

    # Calculate ratio based on selected move
    l_up = flag == 1 ? mc.kappa_up[i] : mc.kappa_up[site]
    l_down = flag == 1 ? mc.kappa_down[site] : mc.kappa_down[i]

    # Calculate acceptance ratio
    ratio = if flag == 1
        mc.W_up[site, l_up] * mc.W_down[i, l_down]
    else
        mc.W_up[i, l_up] * mc.W_down[site, l_down]
    end
    let acc_prob = abs2(ratio)
        if acc_prob >= 1 && r < Zmu / Zmax
            update_configurations!(mc, flag, i, site, l_up, l_down)
            measure!(ctx, :acc, 1.0)
        elseif acc_prob < 1 && r < (Zmu / Zmax) * acc_prob
            update_configurations!(mc, flag, i, site, l_up, l_down)
            measure!(ctx, :acc, 1.0)
        else
            measure!(ctx, :acc, 0.0)
        end
    end

    # Re-evaluate W matrices periodically
    n_occupied = min(count(!iszero, mc.kappa_up), count(!iszero, mc.kappa_down))
    if ctx.sweeps % n_occupied == 0
        try
            reevaluateW!(mc)
        catch e
            if e isa LinearAlgebra.SingularException
                @warn "lu factorization failed, aborting re-evaluation..."
                throw(e)
            end
        end
    end

    return nothing
end

@inline function Carlo.measure!(mc::MC, ctx::MCContext)
    n_occupied = min(count(!iszero, mc.kappa_up), count(!iszero, mc.kappa_down))
    if ctx.sweeps % n_occupied == 0
        OL = getOL(mc, mc.kappa_up, mc.kappa_down)
        measure!(ctx, :OL, OL)
    end
end

"""    
    Carlo.register_evaluables(::Type{MC}, eval::Evaluator, params::AbstractDict)

Register postprocessing evaluators for final analysis after Monte Carlo simulation.

**IMPORTANT: This is purely for postprocessing** - these evaluators are executed 
at the end of the simulation or during merge operations to compute final observables 
from the collected raw measurements. They do NOT affect the Monte Carlo dynamics.

This function defines how to compute physical observables from the raw measured 
data (:OL values) that was collected during the simulation via Carlo.measure!().

# Arguments
- `::Type{MC}`: Monte Carlo type dispatch
- `eval::Evaluator`: Carlo.jl evaluator for postprocessing computations
- `params::AbstractDict`: Simulation parameters containing lattice dimensions

# Registered Evaluables
- `:energy`: Computes the energy per site from local energy measurements
  - Formula: energy = ⟨OL⟩ / ns, where ns = total number of sites
  - Input: Raw :OL measurements collected during simulation
  - Output: Normalized energy per site for final results

# Usage in Carlo.jl Workflow
1. **During simulation**: Carlo.measure!() collects raw :OL data
2. **After simulation**: This function defines how to postprocess the data
3. **Final analysis**: Carlo.jl automatically applies these evaluators to compute final results
4. **Merge operations**: Used when combining data from multiple simulation runs

# Physics Notes
- OL represents the local energy estimator ⟨x|H|ψ_G⟩/⟨x|ψ_G⟩
- Normalization by ns gives energy per site (intensive quantity)
- Essential for comparing results across different system sizes

# Performance Notes
- Called only once during setup, not during simulation loops
- Actual evaluation happens in postprocessing phase
- Minimal computational overhead during Monte Carlo sampling
"""
@inline function Carlo.register_evaluables(
    ::Type{MC},
    eval::Evaluator,
    params::AbstractDict,
)
    n1 = params[:n1]
    n2 = params[:n2]
    ns = n1 * n2 * 3
    evaluate!(eval, :energy, (:OL,)) do OL
        return OL / ns
    end
    return nothing
end

"""    
    Carlo.write_checkpoint(mc::MC, out::HDF5.Group)

Save Monte Carlo state for resuming simulations or postprocessing analysis.

**Purpose**: Serializes the current spinon configuration to an HDF5 file for later 
use in simulation continuation or postprocessing workflows.

# Arguments
- `mc::MC`: Monte Carlo state to save
- `out::HDF5.Group`: Output HDF5 group for writing checkpoint data

# Saved Data
- `kappa_up`: Current up-spinon configuration vector
- `kappa_down`: Current down-spinon configuration vector

# Usage Contexts
1. **Simulation checkpointing**: Save state for restarting long simulations
2. **Postprocessing**: Preserve final configurations for additional analysis
3. **Data archival**: Store representative configurations for later study

# Notes
- Called automatically by Carlo.jl during simulation checkpointing
- Essential for resuming interrupted simulations
- Enables analysis of specific configurations post-simulation
"""
@inline function Carlo.write_checkpoint(mc::MC, out::HDF5.Group)
    out["kappa_up"] = Vector{Int}(mc.kappa_up)
    out["kappa_down"] = Vector{Int}(mc.kappa_down)
    return nothing
end

"""    
    Carlo.read_checkpoint!(mc::MC, in::HDF5.Group)

Restore Monte Carlo state from saved checkpoint for simulation continuation.

**Purpose**: Loads previously saved spinon configurations to resume a simulation 
from a specific state or initialize postprocessing analysis.

# Arguments
- `mc::MC`: Monte Carlo object to restore (modified in-place)
- `in::HDF5.Group`: Input HDF5 group containing checkpoint data

# Restored Data
- `kappa_up`: Up-spinon configuration vector
- `kappa_down`: Down-spinon configuration vector

# Usage Contexts
1. **Simulation resumption**: Continue interrupted long simulations
2. **Postprocessing initialization**: Start analysis from specific configurations
3. **Reproducibility**: Restore exact simulation states

# Side Effects
- Modifies mc.kappa_up and mc.kappa_down in-place
- Overwrites current Monte Carlo configuration completely

# Notes
- Called automatically by Carlo.jl when resuming from checkpoints  
- Essential for maintaining simulation continuity
- Enables reproducible analysis workflows
"""
@inline function Carlo.read_checkpoint!(mc::MC, in::HDF5.Group)
    mc.kappa_up = Vector{Int}(read(in, "kappa_up"))
    mc.kappa_down = Vector{Int}(read(in, "kappa_down"))
    return nothing
end

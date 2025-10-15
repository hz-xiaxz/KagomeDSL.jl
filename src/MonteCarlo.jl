mutable struct MCState{N_up,N_down} <: Carlo.AbstractMC
    Ham::Hamiltonian{N_up,N_down}
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

const MC = MCState

"""
    MCState

Represents the current state of a Variational Monte Carlo (VMC) simulation.
This mutable struct holds all necessary information to perform Monte Carlo
updates and measurements for a given quantum spin liquid system.

# Fields
- `Ham::Hamiltonian{N_up, N_down}`: The Hamiltonian defining the physical system,
  parameterized by the number of up (`N_up`) and down (`N_down`) spinons.
- `kappa_up::Vector{Int}`: The configuration vector for up-spinons. Non-zero
  entries indicate occupied sites, with the value representing the orbital index.
- `kappa_down::Vector{Int}`: The configuration vector for down-spinons, similar to `kappa_up`.
- `W_up::AbstractMatrix`: The one-particle Green's function matrix for up-spinons.
  This matrix is crucial for efficient calculation of wavefunction ratios during updates.
- `W_down::AbstractMatrix`: The one-particle Green's function matrix for down-spinons.
- `W_up_col_cache::AbstractVector`: Pre-allocated cache for column operations during `W_up` updates.
  Used to minimize memory allocations and improve performance.
- `W_up_row_cache::AbstractVector`: Pre-allocated cache for row operations during `W_up` updates.
- `W_down_col_cache::AbstractVector`: Pre-allocated cache for column operations during `W_down` updates.
- `W_down_row_cache::AbstractVector`: Pre-allocated cache for row operations during `W_down` updates.

# Notes
- The `kappa` vectors encode the real-space configuration of spinons.
- The `W` matrices are derived from the Hamiltonian's orbitals and the `kappa` configurations.
- The cache arrays are essential for achieving high performance through in-place rank-1 updates
  of the Green's function matrices, avoiding costly memory reallocations.
"""

"""
    reevaluateW!(mc::MCState)

Recalculate the one-particle Green's function matrices W_up and W_down.

This function reconstructs the Green's functions from the current spin configuration
by solving the linear systems involving the tilde_U matrices.

# Arguments
- `mc::MCState`: Monte Carlo state object

# Algorithm
1. Construct tilde_U matrices from current kappa configurations
2. Solve linear systems: ``tilde_U_up \\ I`` and ``tilde_U_down \\ I``
3. Compute Green's functions: ``W = U * (tilde_U \\ I)``

# Note
``U \\ I`` is numerically more stable than computing the inverse directly.

The Green's function matrices W represent the one-particle propagators
``⟨x|cᵢ⁺cⱼ|ψ⟩/⟨x|ψ⟩`` for the current configuration.

Periodic re-evaluation is necessary to override numerical instability that
accumulates from repeated rank-1 updates to the Green's function matrices.
"""
function reevaluateW!(mc::MCState)
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

Construct the tilde_U matrix by selecting and rearranging rows from U.

The tilde_U matrix is formed by selecting the rows of U corresponding to occupied
sites (where kappa[l] ≠ 0) and placing them in the order specified by kappa.

# Arguments
- `U::AbstractMatrix`: Source matrix of size (n × m) where n is number of sites
- `kappa::Vector{Int}`: Configuration vector where non-zero values indicate occupied sites

# Returns
- A matrix of size (m × m) with same element type as U

# Throws
- `DimensionMismatch`: If length(kappa) ≠ number of rows in U
- `ArgumentError`: If kappa does not contain exactly m non-zero entries
- `BoundsError`: If any non-zero kappa value is out of bounds

# Note
This matrix is used in the construction of the one-particle Green's function W.
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

function relabel_configuration!(kappa::Vector{Int})
    occupied_sites = findall(!iszero, kappa)
    new_kappa = zeros(Int, length(kappa))
    for (new_label, site) in enumerate(occupied_sites)
        new_kappa[site] = new_label
    end
    copyto!(kappa, new_kappa)
    return nothing
end

function spin_plus_transition(mc_n::MCState{N_up,N_down}, site::Int) where {N_up,N_down}
    # 1. Check for a down spin at the site
    if !is_occupied(mc_n.kappa_down, site)
        throw(ArgumentError("No down spin to flip at site $site"))
    end

    # 2. Create new kappa vectors
    new_kappa_up = copy(mc_n.kappa_up)
    new_kappa_down = copy(mc_n.kappa_down)

    # 3. Apply spin flip: site 'site' loses a down spin and gains an up spin
    # The new up spin gets the (N_up + 1)-th orbital
    new_kappa_up[site] = N_up + 1 # Assign the new orbital index
    new_kappa_down[site] = 0

    # 4. Relabel
    relabel_configuration!(new_kappa_up)
    relabel_configuration!(new_kappa_down)

    # 5. Create new Hamiltonian and MCState
    new_N_up = N_up + 1
    new_N_down = N_down - 1

    new_ham = Hamiltonian{new_N_up,new_N_down}(
        mc_n.Ham.U_up_plus,
        mc_n.Ham.U_down_minus,
        mc_n.Ham.H_mat,
        mc_n.Ham.nn,
        # For the new U_up_plus and U_down_minus, we would need orbitals for
        # (N_up+2, N_down-2) and (N_up, N_down) sectors respectively.
        # For now, we'll leave them empty since they're not needed immediately.
        zeros(ComplexF64, size(mc_n.Ham.H_mat, 1), new_N_up + 1),
        zeros(ComplexF64, size(mc_n.Ham.H_mat, 1), max(0, new_N_down - 1)),
        zeros(ComplexF64, size(mc_n.Ham.H_mat, 1), max(0, new_N_up - 1)),
        zeros(ComplexF64, size(mc_n.Ham.H_mat, 1), new_N_down + 1),
    )

    # Create the new MCState. W matrices are zero-initialized.
    new_mc = MCState(
        new_ham,
        new_kappa_up,
        new_kappa_down,
        zeros(eltype(new_ham.U_up), size(new_ham.U_up, 1), new_N_up),
        zeros(eltype(new_ham.U_down), size(new_ham.U_down, 1), new_N_down),
    )

    # The W matrices need to be calculated.
    reevaluateW!(new_mc)

    return new_mc
end

function apply_operator(
    op::SpinPlusOperator{N_up,N_down},
    state::MCState{N_up,N_down},
) where {N_up,N_down}
    return spin_plus_transition(state, op.site)
end

function spin_minus_transition(mc_n::MCState{N_up,N_down}, site::Int) where {N_up,N_down}
    # 1. Check for an up spin at the site
    if !is_occupied(mc_n.kappa_up, site)
        throw(ArgumentError("No up spin to flip at site $site"))
    end

    # 2. Create new kappa vectors
    new_kappa_up = copy(mc_n.kappa_up)
    new_kappa_down = copy(mc_n.kappa_down)

    # 3. Apply spin flip: site 'site' loses an up spin and gains a down spin
    # The new down spin gets the (N_down + 1)-th orbital
    new_kappa_down[site] = N_down + 1 # Assign the new orbital index
    new_kappa_up[site] = 0

    # 4. Relabel
    relabel_configuration!(new_kappa_up)
    relabel_configuration!(new_kappa_down)

    # 5. Create new Hamiltonian and MCState
    new_N_up = N_up - 1
    new_N_down = N_down + 1

    new_ham = Hamiltonian{new_N_up,new_N_down}(
        mc_n.Ham.U_up_minus,
        mc_n.Ham.U_down_plus,
        mc_n.Ham.H_mat,
        mc_n.Ham.nn,
        # For the new U_up_plus and U_down_minus, we would need orbitals for
        # (N_up, N_down) and (N_up-2, N_down+2) sectors respectively.
        # For now, we'll leave them empty since they're not needed immediately.
        zeros(ComplexF64, size(mc_n.Ham.H_mat, 1), new_N_up + 1),
        zeros(ComplexF64, size(mc_n.Ham.H_mat, 1), max(0, new_N_down - 1)),
        zeros(ComplexF64, size(mc_n.Ham.H_mat, 1), max(0, new_N_up - 1)),
        zeros(ComplexF64, size(mc_n.Ham.H_mat, 1), new_N_down + 1),
    )

    # Create the new MCState. W matrices are zero-initialized.
    new_mc = MCState(
        new_ham,
        new_kappa_up,
        new_kappa_down,
        zeros(eltype(new_ham.U_up), size(new_ham.U_up, 1), new_N_up),
        zeros(eltype(new_ham.U_down), size(new_ham.U_down, 1), new_N_down),
    )

    # The W matrices need to be calculated.
    reevaluateW!(new_mc)

    return new_mc
end

function apply_operator(
    op::SpinMinusOperator{N_up,N_down},
    state::MCState{N_up,N_down},
) where {N_up,N_down}
    return spin_minus_transition(state, op.site)
end

function get_log_det_ratio(
    mc_n::MCState{N_up,N_down},
    mc_np1::MCState{N_up_p1,N_down_m1},
) where {N_up,N_down,N_up_p1,N_down_m1}
    # This is inefficient as it re-calculates tilde_U matrices

    # logdet for state n
    tilde_U_up_n = tilde_U(mc_n.Ham.U_up, mc_n.kappa_up)
    tilde_U_down_n = tilde_U(mc_n.Ham.U_down, mc_n.kappa_down)
    log_det_n = real(logdet(tilde_U_up_n)) + real(logdet(tilde_U_down_n))

    # logdet for state n+1
    tilde_U_up_np1 = tilde_U(mc_np1.Ham.U_up, mc_np1.kappa_up)
    tilde_U_down_np1 = tilde_U(mc_np1.Ham.U_down, mc_np1.kappa_down)
    log_det_np1 = real(logdet(tilde_U_up_np1)) + real(logdet(tilde_U_down_np1))

    return log_det_np1 - log_det_n
end

function measure_S_plus(mc::MCState, site::Int)
    if !is_occupied(mc.kappa_down, site)
        return 0.0 + 0.0im, 0.0 # return amplitude and amplitude squared
    end

    mc_np1 = spin_plus_transition(mc, site)

    log_ratio = get_log_det_ratio(mc, mc_np1)

    ratio = complex(exp(log_ratio))
    return ratio, abs2(ratio)
end


"""
    is_occupied(kappa::Vector{Int}, l::Int) -> Bool

Check if site `l` is occupied in the configuration vector.

# Arguments
- `kappa::Vector{Int}`: Configuration vector where non-zero values indicate occupation
- `l::Int`: Site index to check

# Returns
- `Bool`: true if site l is occupied (kappa[l] ≠ 0), false otherwise

# Throws
- `BoundsError`: if l is outside the valid range [1, length(kappa)]
"""
@inline function is_occupied(kappa::Vector{Int}, l::Int)
    @boundscheck 1 ≤ l ≤ length(kappa) || throw(BoundsError(kappa, l))
    @inbounds return kappa[l] != 0
end

"""
    MCState(params::AbstractDict)

Create a Monte Carlo object from a dictionary of parameters.

This is the user-facing, high-level constructor that initializes the Monte Carlo
simulation state with appropriate dimensions and default configurations.

# Arguments
- `params::AbstractDict`: Dictionary containing simulation parameters

# Required Parameters
- `:n1::Int`: Number of unit cells in x-direction
- `:n2::Int`: Number of unit cells in y-direction
- `:PBC::Tuple{Bool,2}`: Periodic boundary conditions
- `:N_up::Int`: Number of up-spin electrons
- `:N_down::Int`: Number of down-spin electrons

# Optional Parameters
- `:antiPBC::Tuple{Bool,2}`: Anti-periodic boundary conditions (default: (false, false))
- `:lattice::Type`: Lattice type (default: DoubleKagome)
- `:B::Float64`: Magnetic field strength (default: 0.0)
- `:link_in::Function`: Intra-cell linking function (default: pi_link_in)
- `:link_inter::Function`: Inter-cell linking function (default: pi_link_inter)

# Returns
- `MCState`: Initialized Monte Carlo state object
"""
function MCState(params::AbstractDict)
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

    return MCState(Ham, kappa_up, kappa_down, W_up, W_down)
end

"""
    MCState(Ham, kappa_up, kappa_down, W_up, W_down)

Create a Monte Carlo object from its core components.

This constructor automatically creates the necessary cache arrays for efficient
matrix updates. It's primarily used for internal logic and testing purposes.

# Arguments
- `Ham::Hamiltonian`: The Hamiltonian describing the physical system
- `kappa_up::Vector{Int}`: Initial up-spin configuration vector
- `kappa_down::Vector{Int}`: Initial down-spin configuration vector
- `W_up::AbstractMatrix`: Initial one-particle Green's function for up-spin
- `W_down::AbstractMatrix`: Initial one-particle Green's function for down-spin

# Returns
- `MCState`: Monte Carlo state object with initialized cache arrays

# Note
The cache arrays enable zero-allocation updates of the Green's function matrices
using rank-1 update operations.
"""
function MCState(
    Ham::Hamiltonian{N_up,N_down},
    kappa_up::Vector{Int},
    kappa_down::Vector{Int},
    W_up::AbstractMatrix,
    W_down::AbstractMatrix,
) where {N_up,N_down}
    ns, N_up_W = size(W_up)
    _, N_down_W = size(W_down)
    @assert ns != 0
    @assert N_up_W != 0
    @assert N_down_W != 0
    @assert N_up_W == N_up
    @assert N_down_W == N_down


    return MCState{N_up,N_down}(
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
    update_W_matrices!(mc::MCState; K_up::Int, K_down::Int, l_up::Int, l_down::Int)

Update the one-particle Green's function matrices using rank-1 updates.

This function performs simultaneous updates of both W_up and W_down matrices
using the Sherman-Morrison formula for efficient matrix inversion updates.

# Arguments
- `mc::MCState`: Monte Carlo state object
- `K_up::Int`: Row index for up-spin matrix update
- `K_down::Int`: Row index for down-spin matrix update
- `l_up::Int`: Column index for up-spin matrix update
- `l_down::Int`: Column index for down-spin matrix update

# Algorithm
Performs rank-1 updates: `W' = W - (W[:,l] * (W[K,:] - δ_{l,:})') / W[K,l]`
for both spin components using pre-allocated cache arrays.
"""
function update_W_matrices!(mc::MCState; K_up::Int, K_down::Int, l_up::Int, l_down::Int)
    update_W!(mc.W_up, l_up, K_up, mc.W_up_col_cache, mc.W_up_row_cache)
    update_W!(mc.W_down, l_down, K_down, mc.W_down_col_cache, mc.W_down_row_cache)
end

"""
    update_W!(W::AbstractMatrix, l::Int, K::Int, col_cache::AbstractVector, row_cache::AbstractVector)

Perform a rank-1 update on the Green's function matrix using Sherman-Morrison formula.

# Arguments
- `W::AbstractMatrix`: Green's function matrix to update
- `l::Int`: Column index for the update
- `K::Int`: Row index for the update
- `col_cache::AbstractVector`: Pre-allocated cache for column operations
- `row_cache::AbstractVector`: Pre-allocated cache for row operations

# Algorithm
Implements the Sherman-Morrison update:
``W'_{I,j} = W_{I,j} - W_{I,l} / W_{K,l} * (W_{K,j} - δ_{l,j})``

This represents a rank-1 update to the matrix inverse that maintains the
Green's function relationship after a single-particle move.
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
    init_conf_qr!(mc::MCState, ns::Int, N_up::Int)

Initialize particle configurations using QR decomposition with column pivoting.

Reference: Quantum Monte Carlo Approaches for Correlated Systems (Becca and Sorella, 2017) P130

As system size increases, ``⟨Φ|x⟩`` becomes exponentially small compared to system size,
inversely proportional to the dimension of the Hilbert space.

The smallness of initial ``⟨Φ|x⟩`` leads to numerical instability in the first computation
of the W matrices, but does not affect Markov chain sampling since we always calculate
``⟨Φ|x'⟩/⟨Φ|x⟩``, which should not be numerically unstable.

This method selects sites that maximize linear independence from the U matrices to
ensure non-singular tilde_U matrices and avoid random trial-and-error initialization.

# Arguments
- `mc::MCState`: Monte Carlo state object (modified in-place)
- `ns::Int`: Total number of sites in the lattice
- `N_up::Int`: Number of up-spin electrons

# Algorithm
1. Perform QR decomposition on `mc.Ham.U_up'` to select `N_up` most linearly independent
   rows (sites) for up-spin electrons
2. Use remaining unoccupied sites and perform QR decomposition on corresponding subset
   of `mc.Ham.U_down'` to select `N_down` sites for down-spin electrons

# Note
This deterministic approach produces non-singular tilde_U matrices for stable Green's function initialization.
"""
function init_conf_qr!(mc::MCState, ns::Int, N_up::Int, N_down::Int)
    # Use QR decomposition with pivoting to select linearly independent rows from U
    # This should give a non-singular tilde_U matrix if possible.

    # For kappa_up
    F_up = qr(mc.Ham.U_up', ColumnNorm())
    sites_up = F_up.p[1:N_up]
    mc.kappa_up = zeros(Int, ns)
    for (i, site) in enumerate(sites_up)
        mc.kappa_up[site] = i
    end

    # For kappa_down
    if N_down > 0
        available_sites = setdiff(1:ns, sites_up)
        U_down_subset = mc.Ham.U_down[collect(available_sites), :]

        F_down = qr(U_down_subset', ColumnNorm())
        sites_down_indices = F_down.p[1:N_down]
        sites_down = collect(available_sites)[sites_down_indices]

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
    find_initial_configuration!(mc::MCState, ns::Int, N_up::Int)

Find a non-singular initial configuration for Monte Carlo simulation.

This function uses a deterministic QR-based method to initialize particle
configurations and compute the initial Green's function matrices.

# Arguments
- `mc::MCState`: Monte Carlo state object (modified in-place)
- `ns::Int`: Total number of sites in the lattice
- `N_up::Int`: Number of up-spin electrons

# Algorithm
1. Initialize configurations using QR decomposition with column pivoting
2. Construct tilde_U matrices from the initialized configurations
3. Solve linear systems to compute Green's function matrices
4. Handle singular matrix exceptions with informative error messages

# Throws
- `ErrorException`: If QR-based configuration produces singular matrices,
  indicating potential rank deficiency in the Hamiltonian
"""
function find_initial_configuration!(mc::MCState, ns::Int, N_up::Int, N_down::Int)
    init_conf_qr!(mc, ns, N_up, N_down)

    tilde_U_up = tilde_U(mc.Ham.U_up, mc.kappa_up)
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
    Carlo.init!(mc::MCState, ctx::MCContext, params::AbstractDict)

Initialize the Monte Carlo object for Carlo.jl framework integration.

This function sets up the Monte Carlo simulation by finding an initial
non-singular configuration using QR-based initialization.

# Arguments
- `mc::MCState`: Monte Carlo state object
- `ctx::MCContext`: Carlo.jl context object
- `params::AbstractDict`: Simulation parameters dictionary

# Parameters Used
- `:n1::Int`: Number of unit cells in x-direction
- `:n2::Int`: Number of unit cells in y-direction
- `:N_up::Int`: Number of up-spin electrons

# Note
This function integrates with the Carlo.jl framework and is called automatically
during simulation initialization to prepare the Monte Carlo state.
"""
@inline function Carlo.init!(mc::MCState, ctx::MCContext, params::AbstractDict)
    n1 = params[:n1]
    n2 = params[:n2]
    ns = n1 * n2 * 3
    N_up = params[:N_up]
    N_down = params[:N_down]
    find_initial_configuration!(mc, ns, N_up, N_down)
end

"""
    Z(nn::AbstractArray, kappa_up::AbstractVector, kappa_down::AbstractVector) -> Int

Count the number of valid spin exchange moves between neighboring sites.

This function calculates the number of neighboring site pairs where one site
is occupied by an up-spin electron and the other by a down-spin electron,
indicating potential spin exchange moves.

# Arguments
- `nn::AbstractArray`: Array of neighbor bonds (site pairs)
- `kappa_up::AbstractVector`: Up-spin configuration vector
- `kappa_down::AbstractVector`: Down-spin configuration vector

# Returns
- `Int`: Number of valid spin exchange move opportunities
"""
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

"""
    update_configurations!(mc, flag::Int, i::Int, site::Int, l_up::Int, l_down::Int)

Update both Green's function matrices and configuration vectors for a spin exchange move.

This function performs the complete update for a spin exchange between sites i and site,
including both the Green's function matrices (W_up, W_down) and the configuration
vectors (kappa_up, kappa_down).

# Arguments
- `mc`: Monte Carlo state object
- `flag::Int`: Move type indicator (1 or 2)
- `i::Int`: First site index
- `site::Int`: Second site index
- `l_up::Int`: Up-spin label to move
- `l_down::Int`: Down-spin label to move

# Move Types
- `flag = 1`: Move up-spin from site i to site, down-spin from site to i
- `flag = 2`: Move up-spin from site to i, down-spin from site i to site
"""
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

    # --- ADDED: Ensure kappa vectors are always canonical ---
    relabel_configuration!(mc.kappa_up)
    relabel_configuration!(mc.kappa_down)
    # -------------------------------------------------------
end

"""
    Carlo.sweep!(mc::MCState, ctx::MCContext) -> Nothing

Perform one Monte Carlo sweep for Mott state simulation.

This function implements a two-spin swap update where electrons can exchange
positions between different spin states at occupied neighboring sites.

# Arguments
- `mc::MCState`: Monte Carlo state object
- `ctx::MCContext`: Carlo.jl context object

# Algorithm
1. Calculate number of valid spin exchange moves Zμ
2. Use Metropolis acceptance criterion based on move probability
3. Select random neighbor bond for potential exchange
4. Calculate acceptance ratio using Green's function matrix elements
5. Update configurations if move is accepted
6. Periodically re-evaluate W matrices to maintain numerical stability

# Note
The W matrices are re-evaluated periodically (every n_occupied sweeps) to
override numerical instability from repeated rank-1 updates. If re-evaluation
fails due to singular matrices, a warning is issued and simulation continues.
"""
@inline function Carlo.sweep!(mc::MCState, ctx::MCContext)
    # MCMC scheme for Mott state
    # Electrons swap between different spins at occupied sites
    nn = mc.Ham.nn
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

"""
    @inline function Carlo.measure!(mc::MCState, ctx::MCContext)

Measures observables during the simulation and collects data.

This function is called periodically during the simulation to measure
physical observables like local energy and collect them for postprocessing.

# Arguments
- `mc::MCState`: Monte Carlo state object
- `ctx::MCContext`: Carlo.jl context object

# Measured Quantities
- `:OL`: Local energy estimator `getOL(mc, mc.kappa_up, mc.kappa_down)`

# Note
Measurements are taken periodically (every n_occupied sweeps) to reduce
correlation between samples and improve measurement efficiency.
"""
@inline function Carlo.measure!(mc::MCState, ctx::MCContext)
    n_occupied = min(count(!iszero, mc.kappa_up), count(!iszero, mc.kappa_down))
    if ctx.sweeps % n_occupied == 0
        OL = getOL(mc, mc.kappa_up, mc.kappa_down)
        measure!(ctx, :OL, OL)

        # New measurement
        ns = length(mc.kappa_up)
        s_plus_amp = 0.0 + 0.0im
        s_plus_amp_sq = 0.0
        for i = 1:ns
            amp, amp_sq = measure_S_plus(mc, i)
            s_plus_amp += amp
            s_plus_amp_sq += amp_sq
        end
        measure!(ctx, :S_plus_amp, s_plus_amp / ns)
        measure!(ctx, :S_plus_amp_sq, s_plus_amp_sq / ns)
    end
end

"""
    Carlo.register_evaluables(::Type{MCState}, eval::Evaluator, params::AbstractDict)

Register postprocessing evaluators for final analysis after Monte Carlo simulation.

**IMPORTANT: This is purely for postprocessing** - these evaluators are executed
at the end of the simulation or during merge operations to compute final observables
from the collected raw measurements. They do NOT affect the Monte Carlo dynamics.

This function defines how to compute physical observables from the raw measured
data (:OL values) that was collected during the simulation via Carlo.measure!().

# Arguments
- `::Type{MCState}`: Monte Carlo type dispatch
- `eval::Evaluator`: Carlo.jl evaluator for postprocessing computations
- `params::AbstractDict`: Simulation parameters containing lattice dimensions

# Registered Evaluables
- `:energy`: Computes the energy per site from local energy measurements
  - Formula: ``energy = ⟨OL⟩ / ns``, where ns = total number of sites
  - Input: Raw :OL measurements collected during simulation
  - Output: Normalized energy per site for final results

# Usage in Carlo.jl Workflow
1. **During simulation**: `Carlo.measure!()` collects raw `:OL` data
2. **After simulation**: This function defines how to postprocess the data
3. **Final analysis**: Carlo.jl automatically applies these evaluators to compute final results
4. **Merge operations**: Used when combining data from multiple simulation runs

# Physics Notes
- `OL` represents the local energy estimator ``⟨x|H|ψ_G⟩/⟨x|ψ_G⟩``
- Normalization by `ns` gives energy per site (intensive quantity)
- Essential for comparing results across different system sizes

# Performance Notes
- Called only once during setup, not during simulation loops
- Actual evaluation happens in postprocessing phase
- Minimal computational overhead during Monte Carlo sampling
"""
@inline function Carlo.register_evaluables(
    ::Type{<:MCState},
    eval::Evaluator,
    params::AbstractDict,
)
    n1 = params[:n1]
    n2 = params[:n2]
    ns = n1 * n2 * 3
    evaluate!(eval, :energy, (:OL,)) do OL
        OL / ns
    end
    evaluate!(eval, :S_plus_amp, (:S_plus_amp,)) do S_plus_amp
        S_plus_amp
    end
    evaluate!(eval, :S_plus_amp_sq, (:S_plus_amp_sq,)) do S_plus_amp_sq
        S_plus_amp_sq
    end
end

"""
    Carlo.write_checkpoint(mc::MCState, out::HDF5.Group)

Save Monte Carlo state for resuming simulations or postprocessing analysis.

**Purpose**: Serializes the current spin configuration to an HDF5 file for later
use in simulation continuation or postprocessing workflows.

# Arguments
- `mc::MCState`: Monte Carlo state to save
- `out::HDF5.Group`: Output HDF5 group for writing checkpoint data

# Saved Data
- `kappa_up`: Current up-spin configuration vector
- `kappa_down`: Current down-spin configuration vector

# Usage Contexts
1. **Simulation checkpointing**: Save state for restarting long simulations
2. **Postprocessing**: Preserve final configurations for additional analysis
3. **Data archival**: Store representative configurations for later study

# Notes
- Called automatically by Carlo.jl during simulation checkpointing
- Essential for resuming interrupted simulations
- Enables analysis of specific configurations post-simulation
"""
@inline function Carlo.write_checkpoint(mc::MCState, out::HDF5.Group)
    out["kappa_up"] = mc.kappa_up
    out["kappa_down"] = mc.kappa_down
    return nothing
end

"""
    Carlo.read_checkpoint!(mc::MCState, in::HDF5.Group)

Restore Monte Carlo state from saved checkpoint for simulation continuation.

**Purpose**: Loads previously saved spin configurations to resume a simulation
from a specific state or initialize postprocessing analysis.

# Arguments
- `mc::MCState`: Monte Carlo object to restore (modified in-place)
- `in::HDF5.Group`: Input HDF5 group containing checkpoint data

# Restored Data
- `kappa_up`: Up-spin configuration vector
- `kappa_down`: Down-spin configuration vector

# Usage Contexts
1. **Simulation resumption**: Continue interrupted long simulations
2. **Postprocessing initialization**: Start analysis from specific configurations
3. **Reproducibility**: Restore exact simulation states

# Side Effects
- Modifies `mc.kappa_up` and `mc.kappa_down` in-place
- Overwrites current Monte Carlo configuration completely

# Notes
- Called automatically by Carlo.jl when resuming from checkpoints
- Essential for maintaining simulation continuity
- Enables reproducible analysis workflows
"""
@inline function Carlo.read_checkpoint!(mc::MCState, in::HDF5.Group)
    mc.kappa_up = read(in, "kappa_up")
    mc.kappa_down = read(in, "kappa_down")
    return nothing
end

@doc raw"""
    getOL(mc::AbstractMC, kappa_up, kappa_down) -> Float64

Compute the local energy estimator for quantum Monte Carlo.

Calculates the observable O_L = ⟨x|H|ψ_G⟩/⟨x|ψ_G⟩, which provides an unbiased
estimator of the ground state energy in the Stochastic Series Expansion method.

In the spinon mean-field framework:
- |ψ_G⟩ is the Slater Determinant spinon state (Gutzwiller projected)
- |x⟩ = |κ_up⟩ ⊗ |κ_down⟩ is a particular spinon configuration
- H is the original Heisenberg Hamiltonian (not the mean-field one!)

The ratio gives the local contribution to the energy from configuration |x⟩.

# Arguments
- `mc::AbstractMC`: Monte Carlo state containing W_up, W_down matrices
- `kappa_up::Vector{Int}`: Current up-spinon configuration
- `kappa_down::Vector{Int}`: Current down-spinon configuration

# Returns
- `Float64`: Local energy contribution from this configuration

# Algorithm
1. Compute H|x⟩ using getxprime() to get all reachable configurations
2. For diagonal terms: directly add the Ising contributions
3. For off-diagonal terms: weight by the wavefunction amplitude ratios W_up, W_down
4. Sum all contributions to get the local energy

# Physics Notes
- This is the "local energy" in variational Monte Carlo terminology
- Fluctuations in O_L reflect the quality of the trial wavefunction
- Used to measure energy and other observables in quantum spin liquids
- The Hamiltonian H must be the original physical one, not the mean-field approximation

# Performance Notes
- Critical function called in every Monte Carlo step
- Bounds checking disabled with @inbounds for performance
- Uses efficient iteration over pairs() for dictionary access
"""
@inline function getOL(mc::MCState, kappa_up::Vector{Int}, kappa_down::Vector{Int})
    @assert length(kappa_up) == length(kappa_down) "The length of the up and down configurations should be the same, got: $(length(kappa_up)) and $(length(kappa_down))"
    # if double occupied state, no possibility to have a non-zero overlap

    OL = 0.0
    xprime = getxprime(mc.Ham, kappa_up, kappa_down)
    @inbounds for (conf, coff) in pairs(xprime)
        if conf == (-1, -1, -1, -1)
            OL += coff
        else
            update_up = mc.W_up[conf[1], conf[2]]
            update_down = mc.W_down[conf[3], conf[4]]
            OL += coff * update_up * update_down
        end
    end
    return OL
end

mutable struct MC <: AbstractMC
    Ham::Hamiltonian
    kappa_up::Vector{Int}
    kappa_down::Vector{Int}
    W_up::AbstractMatrix
    W_down::AbstractMatrix
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

function init_conf(rng::AbstractRNG, ns::Int, N_up::Int)
    # dealing with conf_up
    # Initialize array with zeros
    kappa_up = zeros(Int, ns)

    # Generate random positions for 1:N_up
    positions = randperm(rng, ns)[1:N_up]

    # Place 1:N_up at the random positions
    for (i, pos) in enumerate(positions)
        kappa_up[pos] = i
    end
    kappa_down = zeros(Int, ns)
    # kappa_down should occupy the rest of the sites
    res_pos = randperm(rng, ns - N_up)
    for i = 1:ns
        if kappa_up[i] == 0
            kappa_down[i] = pop!(res_pos)
        end
    end

    return kappa_up, kappa_down
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
    @inbounds return !iszero(kappa[l])
end

"""
    MC(params::AbstractDict)
------------
Create a Monte Carlo object
"""
function MC(params::AbstractDict)
    n1 = params[:n1]
    n2 = params[:n2]
    PBC = params[:PBC]
    lat = DoubleKagome(1.0, n1, n2, PBC)
    N_up = params[:N_up]
    N_down = params[:N_down]
    Ham = Hamiltonian(N_up, N_down, lat)
    rng = Random.Xoshiro(42)
    ns = n1 * n2 * 3
    kappa_up, kappa_down = init_conf(rng, ns, N_up)
    tilde_U_up = tilde_U(Ham.U_up, kappa_up)
    tilde_U_down = tilde_U(Ham.U_down, kappa_down)
    U_upinvs = tilde_U_up \ I
    U_downinvs = tilde_U_down \ I
    # calculate W_up and W_down
    # Calculate W matrices using matrix multiplication
    W_up = Ham.U_up * U_upinvs
    W_down = Ham.U_down * U_downinvs
    return MC(Ham, kappa_up, kappa_down, W_up, W_down)
end

"""
    update_W_matrices(mc::MC; K_up::Int, K_down::Int, l_up::Int, l_down::Int)
------------
Update the W matrices
"""
function update_W_matrices(mc::MC; K_up::Int, K_down::Int, l_up::Int, l_down::Int)
    update_W!(mc.W_up; l = l_up, K = K_up)
    update_W!(mc.W_down; l = l_down, K = K_down)
end

"""
    update_W!(W::AbstractMatrix; l::Int, K::Int)
------------
Update the W matrix
``W'_{I,j} = W_{I,j} - W_{I,l} / W_{K,l} * (W_{K,j} - \\delta_{l,j})``
"""
function update_W!(W::AbstractMatrix; l::Int, K::Int)
    factors = W[:, l] ./ W[K, l]
    Krow = copy(W[K, :])
    Krow[l] -= 1.0
    @views for j in axes(W, 2)
        W[:, j] .-= factors .* Krow[j]
    end
    return nothing
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
    kappa_up, kappa_down = init_conf(ctx.rng, ns, N_up)
    # calculate W_up and W_down
    reevaluateW!(mc)
end

function Z(nn::AbstractArray, kappa_up::AbstractVector, kappa_down::AbstractVector)
    # iterate over all possible moves
    # if two sites connected by a bond, check if they are occupied by different spins
    count = 0
    for bond in nn
        site1 = bond[1]
        site2 = bond[2]
        if is_occupied(kappa_up, site1) && is_occupied(kappa_down, site2)
            count += 1
        elseif is_occupied(kappa_up, site2) && is_occupied(kappa_down, site1)
            count += 1
        end
    end
    return count
end

function update_configurations!(mc, flag::Int, i::Int, site::Int, l_up::Int, l_down::Int)
    if flag == 1
        # Update W matrices
        update_W_matrices(mc; K_up = site, K_down = i, l_up = l_up, l_down = l_down)
        # Update kappa configurations
        mc.kappa_up[i], mc.kappa_up[site] = 0, l_up
        mc.kappa_down[i], mc.kappa_down[site] = l_down, 0
    else
        # Update W matrices
        update_W_matrices(mc; K_up = i, K_down = site, l_up = l_up, l_down = l_down)
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
    if ratio^2 >= 1 && r < Zmu / Zmax
        update_configurations!(mc, flag, i, site, l_up, l_down)
        measure!(ctx, :acc, 1.0)
    elseif ratio^2 < 1 && r < (Zmu / Zmax) * ratio^2
        update_configurations!(mc, flag, i, site, l_up, l_down)
        measure!(ctx, :acc, 1.0)
    else
        measure!(ctx, :acc, 0.0)
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
        measure!(ctx, :OL, OL / 2)
    end
end

@inline function Carlo.register_evaluables(
    ::Type{MC},
    eval::Evaluator,
    params::AbstractDict,
)

    return nothing
end

@inline function Carlo.write_checkpoint(mc::MC, out::HDF5.Group)
    out["kappa_up"] = Vector{Int}(mc.kappa_up)
    out["kappa_down"] = Vector{Int}(mc.kappa_down)
    return nothing
end

@inline function Carlo.read_checkpoint!(mc::MC, in::HDF5.Group)
    mc.kappa_up = Vector{Int}(read(in, "kappa_up"))
    mc.kappa_down = Vector{Int}(read(in, "kappa_down"))
    return nothing
end

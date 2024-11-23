using Random

mutable struct MC <: AbstractMC
    Ham::Hamiltonian
    kappa_up::Vector{Int}
    kappa_down::Vector{Int}
    W_up::AbstractMatrix
    W_down::AbstractMatrix
end

function Carlo.is_thermalized(ctx::MCContext)
    # Hack
    return ctx.thermalization_sweeps == -1
end

function reevaluateW!(mc::MC)
    # Calculate inverse matrices
    tiled_U_up = tiled_U(mc.Ham.U_up, mc.kappa_up)
    tiled_U_down = tiled_U(mc.Ham.U_down, mc.kappa_down)
    let U_upinvs = nothing, U_downinvs = nothing
        try
            U_upinvs = tiled_U_up \ I
        catch e
            if e isa LinearAlgebra.SingularException
                @warn "lu factorization failed, aborting re-evaluation..."
                @show tiled_U_up
            end
        end
        try
            U_downinvs = tiled_U_down \ I
        catch e
            if e isa LinearAlgebra.SingularException
                @warn "lu factorization failed, aborting re-evaluation..."
                @show tiled_U_down
            end
        end

        mc.W_up = mc.Ham.U_up * U_upinvs
        mc.W_down = mc.Ham.U_down * U_downinvs
    end
    # Calculate W matrices using matrix multiplication

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
    tiled_U(U::AbstractMatrix, kappa::Vector{Int})
------------
Creates a tiled matrix by rearranging rows of U according to kappa indices.

Parameters:
- `U`: Source matrix of size (n × m)
- `kappa`: Vector of indices where each non-zero value l indicates that row Rl of U
          should be placed at row l of the output

Returns:
- A matrix of size (m × m) with same element type as U
"""
function tiled_U(U::AbstractMatrix, kappa::Vector{Int})
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
    tiled_U = zeros(eltype(U), m, m)

    @inbounds for (Rl, l) in enumerate(kappa)
        if l != 0
            (1 ≤ l ≤ m) || throw(BoundsError(tiled_U, (l, :)))
            tiled_U[l, :] = U[Rl, :]
        end
    end

    return tiled_U
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
    ensure_nonzero_det(rng, Ham, ns, N_up; max_attempts=100)

Generate configurations with non-zero determinants for both up and down states.

# Arguments
- `rng`: Random number generator
- `Ham`: Hamiltonian object
- `ns`: Number of sites
- `N_up`: Number of up spins
- `max_attempts`: Maximum number of initialization attempts

# Returns
- `kappa_up`, `kappa_down`: Valid configurations with non-zero determinants
"""
function ensure_nonzero_det(rng, Ham, ns, N_up; max_attempts = 100)
    attempts = 0
    kappa_up, kappa_down = init_conf(rng, ns, N_up)

    while attempts < max_attempts
        tiled_U_up = tiled_U(Ham.U_up, kappa_up)
        tiled_U_down = tiled_U(Ham.U_down, kappa_down)

        if abs(det(tiled_U_up)) > 1e-14 && abs(det(tiled_U_down)) > 1e-14
            return kappa_up, kappa_down
        end

        # @warn "det(tiled_U_up) = $(det(tiled_U_up)), det(tiled_U_down) = $(det(tiled_U_down))"
        kappa_up, kappa_down = init_conf(rng, ns, N_up)
        attempts += 1
    end

    error("Failed to find valid configuration after $max_attempts attempts")
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
    kappa_up, kappa_down = ensure_nonzero_det(rng, Ham, ns, N_up, max_attempts = 10000000)
    U_upinvs = tiled_U(Ham.U_up, kappa_up) \ I
    U_downinvs = tiled_U(Ham.U_down, kappa_down) \ I
    # calculate W_up and W_down
    # Calculate W matrices using matrix multiplication
    W_up = Ham.U_up * U_upinvs
    W_down = Ham.U_down * U_downinvs
    return MC(Ham, kappa_up, kappa_down, W_up, W_down)
end

"""
    update_W(W::AbstractMatrix; l::Int, K::Int)
------------
Update the W matrix
``W'_{I,j} = W_{I,j} - W_{I,l} / W_{K,l} * (W_{K,j} - \\delta_{l,j})``
"""
function update_W(W::AbstractMatrix; l::Int, K::Int)
    factors = [W[I, l] / W[K, l] for I in axes(W, 1)]
    new_W = similar(W)
    @inbounds for I in axes(W, 1)
        for j in axes(W, 2)
            new_W[I, j] = W[I, j] - factors[I] * (W[K, j] - (l == j ? 1.0 : 0.0))
        end
    end
    return new_W
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
    ctx.thermalization_sweeps = 0
    mc.kappa_up, mc.kappa_down = ensure_nonzero_det(ctx.rng, mc.Ham, ns, N_up)
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
        mc.W_up = update_W(mc.W_up; l = l_up, K = site)
        mc.W_down = update_W(mc.W_down; l = l_down, K = i)
        # Update kappa configurations
        mc.kappa_up[i], mc.kappa_up[site] = 0, l_up
        mc.kappa_down[i], mc.kappa_down[site] = l_down, 0
    else
        # Update W matrices
        mc.W_up = update_W(mc.W_up; l = l_up, K = i)
        mc.W_down = update_W(mc.W_down; l = l_down, K = site)
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
    if !is_thermalized(ctx)
        det_up = det(tiled_U(mc.Ham.U_up, mc.kappa_up))
        det_down = det(tiled_U(mc.Ham.U_down, mc.kappa_down))
        if abs(det_up) > 1e-14 || abs(det_down) > 1e-14
            ctx.thermalization_sweeps = -1
            @show "thermalized"
        end
    end

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
    # Metropolis acceptance step
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
        let OL = nothing
            function try_measure()
                OL = getOL(mc, mc.kappa_up, mc.kappa_down)
                measure!(ctx, :OL, OL)
                return true
            end

            # First attempt
            try
                return try_measure()
            catch e
                if e isa LinearAlgebra.SingularException
                    @warn "First lu factorization failed, attempting recovery..."
                    # Recovery attempt
                    reevaluateW!(mc)
                    try
                        return try_measure()
                    catch e2
                        if e2 isa LinearAlgebra.SingularException
                            @warn "Recovery failed: lu factorization failed twice"
                            return false
                        end
                        rethrow(e2)  # Rethrow non-singular exceptions
                    end
                end
                rethrow(e)  # Rethrow non-singular exceptions
            end
        end
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

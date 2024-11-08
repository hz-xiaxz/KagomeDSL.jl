mutable struct MC <: AbstractMC
    Ham::Hamiltonian
    conf_up::BitVector
    conf_down::BitVector
    W_up::AbstractMatrix
    W_down::AbstractMatrix
end

function reevaluateW!(mc::MC)
    # Calculate inverse matrices
    U_up_inv = mc.Ham.U_up[Bool.(mc.conf_up), :] \ I
    U_down_inv = mc.Ham.U_down[Bool.(mc.conf_down), :] \ I

    # Calculate W matrices using matrix multiplication
    mc.W_up = mc.Ham.U_up * U_up_inv
    mc.W_down = mc.Ham.U_down * U_down_inv

    return nothing
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
    init_conf = zeros(Bool, ns)
    init_conf[randperm(rng, ns)[1:N_up]] .= true
    conf_up = BitVector(init_conf)
    conf_down = fill(true, ns) - conf_up
    U_upinvs = Ham.U_up[Bool.(conf_up), :] \ I
    U_downinvs = Ham.U_down[Bool.(conf_down), :] \ I
    # calculate W_up and W_down
    # Calculate W matrices using matrix multiplication
    W_up = Ham.U_up * U_upinvs
    W_down = Ham.U_down * U_downinvs
    return MC(Ham, conf_up, conf_down, W_up, W_down)
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
    mc.conf_up = zeros(Bool, ns)
    mc.conf_up[randperm(ctx.rng, ns)[1:N_up]] .= true
    mc.conf_up = BitVector(mc.conf_up)
    mc.conf_down = fill(true, ns) - mc.conf_up
    # calculate W_up and W_down
    reevaluateW!(mc)
end

function getNeigh(rng, ns::Int, nn::AbstractArray)
    while true
        i = rand(rng, 1:ns)
        neigh = filter(x -> x[1] == i, nn)
        if length(neigh) == 0
            continue
        else
            site = sample(rng, neigh)[2]
            return i, site
        end
    end
end

@inline function Carlo.sweep!(mc::MC, ctx::MCContext)
    # should perform MCMC here
    # MCMC scheme, generate a Mott state
    # the electrons are only allowed to swap between different spins and from an occupied site to an unoccupied site
    nn = mc.Ham.nn
    ns = length(mc.conf_up)

    # randomly select a site
    # flip two spins only to perform fast_update
    i, site = getNeigh(ctx.rng, ns, nn)
    ratio = 0

    update_pool = [mc.conf_up[i] && mc.conf_down[site], mc.conf_up[site] && mc.conf_down[i]]
    maybe_update = findall(x -> x == true, update_pool)
    if isempty(maybe_update)
        measure!(ctx, :acc, 0.0)
        return Nothing
    else
        flag = sample(ctx.rng, maybe_update)
        if flag == 1
            # mc.conf_up[i] = false # the old site is empty
            # mc.conf_up[site] = true # the new site is occupied
            # mc.conf_down[i] = true
            # mc.conf_down[site] = false
            l_up = sum(mc.conf_up[1:i])
            l_down = sum(mc.conf_down[1:site])
            ratio = mc.W_up[site, l_up] * mc.W_down[i, l_down]
        elseif flag == 2
            l_up = sum(mc.conf_up[1:site])
            l_down = sum(mc.conf_down[1:i])
            # mc.conf_up[i] = true
            # mc.conf_up[site] = false
            # mc.conf_down[i] = false
            # mc.conf_down[site] = true
            ratio = mc.W_up[i, l_up] * mc.W_down[site, l_down]
        end
    end
    if ratio^2 < 1.0 && rand(ctx.rng) > ratio^2
        # measure acceptance rate here
        measure!(ctx, :acc, 0.0)
    else
        if flag == 1
            mc.W_up = update_W(mc.W_up; l = l_up, K = site)
            mc.W_down = update_W(mc.W_down; l = l_down, K = i)
            mc.conf_up[i] = false
            mc.conf_up[site] = true
            mc.conf_down[i] = true
            mc.conf_down[site] = false
        elseif flag == 2
            mc.W_up = update_W(mc.W_up; l = l_up, K = i)
            mc.W_down = update_W(mc.W_down; l = l_down, K = site)
            mc.conf_up[i] = true
            mc.conf_up[site] = false
            mc.conf_down[i] = false
            mc.conf_down[site] = true
        end
        measure!(ctx, :acc, 1.0)
    end

    # re-evaluate W matrices every Ne sweeps
    number = min(sum(mc.conf_up), sum(mc.conf_down))
    if ctx.sweeps % number == 0
        reevaluateW!(mc)
    end

    return nothing
end

@inline function Carlo.measure!(mc::MC, ctx::MCContext)
    let OL = nothing
        function try_measure()
            OL = getOL(mc, mc.conf_up, mc.conf_down)
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

@inline function Carlo.register_evaluables(
    ::Type{MC},
    eval::Evaluator,
    params::AbstractDict,
)

    return nothing
end

@inline function Carlo.write_checkpoint(mc::MC, out::HDF5.Group)
    out["conf_up"] = Vector{Bool}(mc.conf_up)
    out["conf_down"] = Vector{Bool}(mc.conf_down)
    return nothing
end

@inline function Carlo.read_checkpoint!(mc::MC, in::HDF5.Group)
    mc.conf_up = BitVector(read(in, "conf_up"))
    mc.conf_down = BitVector(read(in, "conf_down"))
    return nothing
end

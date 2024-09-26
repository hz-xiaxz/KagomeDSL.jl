mutable struct MC <: AbstractMC
    Ham::Hamiltonian
    conf_up::BitVector
    conf_down::BitVector
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
    ns = n1 * n2 * 3
    init_conf = zeros(Bool, ns)
    init_conf[randperm(rng, ns)[1:N_up]] .= true
    conf_up = BitVector(init_conf)
    conf_down = fill(true, ns) - conf_up
    return MC(Ham, conf_up, conf_down)
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
    ctx.rng = Random.Xoshiro(42)
    ns = n1 * n2 * 3
    N_up = params[:N_up]
    mc.conf_up = zeros(Bool, ns)
    mc.conf_up[randperm(ctx.rng, ns)[1:N_up]] .= true
    mc.conf_up = BitVector(mc.conf_up)
    mc.conf_down = fill(true, ns) - mc.conf_up
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

function getRatio(
    U::AbstractMatrix,
    Uinvs::AbstractMatrix,
    oldconf::BitVector,
    i::Int,
    site::Int,
)
    l = sum(oldconf[1:i])
    return sum(U[site, :] .* Uinvs[:, l])
end

@inline function Carlo.sweep!(mc::MC, ctx::MCContext)
    # should perform MCMC here
    # MCMC scheme, generate a Mott state
    # the electrons are only allowed to swap between different spins and from an occupied site to an unoccupied site
    nn = mc.Ham.nn
    ns = length(mc.conf_up)
    oldconfup = copy(mc.conf_up)
    oldconfdown = copy(mc.conf_down)
    U_upinvs = mc.Ham.U_up[oldconfup, :] \ I
    U_downinvs = mc.Ham.U_down[oldconfdown, :] \ I

    # randomly select a site
    # flip two spins only to perform fast_update
    i, site = getNeigh(ctx.rng, ns, nn)
    ratio = 0
    if mc.conf_up[i] && mc.conf_down[site]
        mc.conf_up[i] = false # the old site is empty
        mc.conf_up[site] = true # the new site is occupied
        mc.conf_down[i] = true
        mc.conf_down[site] = false
        ratio =
            getRatio(mc.Ham.U_up, U_upinvs, oldconfup, i, site) *
            getRatio(mc.Ham.U_down, U_downinvs, oldconfdown, site, i)
    elseif mc.conf_up[site] && mc.conf_down[i]
        mc.conf_up[i] = true
        mc.conf_up[site] = false
        mc.conf_down[i] = false
        mc.conf_down[site] = true
        ratio =
            getRatio(mc.Ham.U_up, U_upinvs, oldconfup, site, i) *
            getRatio(mc.Ham.U_down, U_downinvs, oldconfdown, i, site)
    elseif mc.conf_up[i] && !mc.conf_up[site]
        # i is occupied, site is not for the up spin
        mc.conf_up[i] = false
        mc.conf_up[site] = true
        ratio = getRatio(mc.Ham.U_up, U_upinvs, oldconfup, i, site)
    elseif !mc.conf_up[i] && mc.conf_up[site]
        # i is empty, site is occupied for the up spin
        mc.conf_up[i] = true
        mc.conf_up[site] = false
        ratio = getRatio(mc.Ham.U_up, U_upinvs, oldconfup, site, i)
    elseif mc.conf_down[i] && !mc.conf_down[site]
        # i is occupied, site is not for the down spin
        mc.conf_down[i] = false
        mc.conf_down[site] = true
        ratio = getRatio(mc.Ham.U_down, U_downinvs, oldconfdown, i, site)
    elseif !mc.conf_down[i] && mc.conf_down[site]
        # i is empty, site is occupied for the down spin
        mc.conf_down[i] = true
        mc.conf_down[site] = false
        ratio = getRatio(mc.Ham.U_down, U_downinvs, oldconfdown, site, i)
    end
    if ratio^2 < 1.0 && rand(ctx.rng) > ratio^2
        mc.conf_up = oldconfup
        mc.conf_down = oldconfdown
    end
    return nothing
end

@inline function Carlo.measure!(mc::MC, ctx::MCContext)
    # get E
    OL = getOL(mc.Ham, mc.conf_up, mc.conf_down)
    measure!(ctx, :OL, OL)
    return nothing
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

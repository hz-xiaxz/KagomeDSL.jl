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
    χ = params[:χ]
    Ham = Hamiltonian(χ, N_up, N_down, lat)
    rng = Random.Xoshiro(42)
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
* `χ` : `Float64` hopping strength
* `N_up` : `Int` number of up spinons
* `N_down` : `Int` number of down spinons
"""
@inline function Carlo.init!(mc::MC, ctx::MCContext, params::AbstractDict)
    n1 = params[:n1]
    n2 = params[:n2]
    ctx.rng = Random.Xoshiro(42)
    ns = n1 * n2 * 3
    conf_up = zeros(Bool, ns)
    N_up = params[:N_up]
    conf_up[randperm(ctx.rng, ns)[1:N_up]] .= true
    mc.conf_up = BitVector(conf_up)
    mc.conf_down = fill(true, ns) - conf_up
end

@inline function Carlo.sweep!(mc::MC, ctx::MCContext)
    # should perform MCMC here
    # MCMC scheme, generate a Mott state
    # the electrons are only allowed to swap between different spins and from an occupied site to an unoccupied site
    nn = mc.Ham.nn
    ns = length(mc.conf_up)
    oldconfup = copy(mc.conf_up)
    oldconfdown = copy(mc.conf_down)

    oldconfupstr = LongBitStr(mc.conf_up)
    oldconfdownstr = LongBitStr(mc.conf_down)
    for i = 1:ns
        neigh = filter(x -> x[1] == i, nn)
        if length(neigh) == 0
            continue
        end
        site = sample(ctx.rng, neigh)[2]
        if mc.conf_up[i] && mc.conf_down[site]
            mc.conf_up[i] = false
            mc.conf_up[site] = true
            mc.conf_down[i] = true
            mc.conf_down[site] = false
        elseif mc.conf_up[site] && mc.conf_down[i]
            mc.conf_up[i] = true
            mc.conf_up[site] = false
            mc.conf_down[i] = false
            mc.conf_down[site] = true
        end
    end
    ratio =
        det(mc.Ham.U_up[mc.conf_up, :]) / det(mc.Ham.U_up[oldconfup, :]) *
        det(mc.Ham.U_down[mc.conf_down, :]) / det(mc.Ham.U_down[oldconfdown, :])
    if ratio < 1.0 || rand(ctx.rng) > ratio
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

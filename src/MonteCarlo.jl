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
    ns = n1 * n2 * 3
    rng = Random.Xoshiro(42)
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
    ns = n1 * n2 * 3
    N_up = params[:N_up]
    N_down = params[:N_down]
    mc.conf_up = zeros(Bool, ns)
    mc.conf_up[randperm(ctx.rng, ns)[1:N_up]] .= true
    mc.conf_up = BitVector(mc.conf_up)
    mc.conf_down = fill(true, ns) - mc.conf_up
end

function getRatio(
    U::AbstractMatrix,
    Uinvs::AbstractMatrix,
    oldconf::BitVector;
    old::Int,
    new::Int,
)
    l = sum(oldconf[1:old])
    new_prime = sum(oldconf[1:new])
    return dot(U[new, :], Uinvs[:, l]) * (-1)^(new_prime - l)
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
    oldconfup = copy(Bool.(mc.conf_up))
    oldconfdown = copy(Bool.(mc.conf_down))
    U_upinvs = mc.Ham.U_up[oldconfup, :] \ I
    U_downinvs = mc.Ham.U_down[oldconfdown, :] \ I

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
            # @assert mc.conf_up[i] && mc.conf_down[site]
            l_up = sum(mc.conf_up[1:i])
            l_down = sum(mc.conf_down[1:site])
            ratio_get_up = getRatio(mc.Ham.U_up, U_upinvs, oldconfup; old = i, new = site)
            new_conf = copy(Bool.(mc.conf_up))
            new_conf[i] = false
            new_conf[site] = true
            ratio_true_up =
                det(mc.Ham.U_up[Bool.(new_conf), :]) / det(mc.Ham.U_up[oldconfup, :])
            if !isapprox(ratio_get_up, ratio_true_up; atol = 1e-14)
                @info ratio_get_up, ratio_true_up
                @show mc.Ham.U_up[site, :] .* U_upinvs[:, l_up]
                @show mc.Ham.U_up[site, :]
                @show U_upinvs[:, l_up]
            end
            ratio =
                getRatio(mc.Ham.U_down, U_downinvs, oldconfdown; old = site, new = i) *
                getRatio(mc.Ham.U_up, U_upinvs, oldconfup; old = i, new = site)
        elseif flag == 2
            # @assert mc.conf_up[site] && mc.conf_down[i]
            l_up = sum(mc.conf_up[1:site])
            l_down = sum(mc.conf_down[1:i])
            # mc.conf_up[i] = true
            # mc.conf_up[site] = false
            # mc.conf_down[i] = false
            # mc.conf_down[site] = true
            ratio =
                getRatio(mc.Ham.U_down, U_downinvs, oldconfdown; old = i, new = site) *
                getRatio(mc.Ham.U_up, U_upinvs, oldconfup, old = site, new = i)
        end
        # elseif mc.conf_up[i] && !mc.conf_up[site]
        #     # i is occupied, site is not for the up spin
        #     mc.conf_up[i] = false
        #     mc.conf_up[site] = true
        #     ratio = getRatio(mc.Ham.U_up, U_upinvs, oldconfup, i, site)
        # elseif !mc.conf_up[i] && mc.conf_up[site]
        #     # i is empty, site is occupied for the up spin
        #     mc.conf_up[i] = true
        #     mc.conf_up[site] = false
        #     ratio = getRatio(mc.Ham.U_up, U_upinvs, oldconfup, site, i)
        # elseif mc.conf_down[i] && !mc.conf_down[site]
        #     # i is occupied, site is not for the down spin
        #     mc.conf_down[i] = false
        #     mc.conf_down[site] = true
        #     ratio = getRatio(mc.Ham.U_down, U_downinvs, oldconfdown, i, site)
        # elseif !mc.conf_down[i] && mc.conf_down[site]
        #     # i is empty, site is occupied for the down spin
        #     mc.conf_down[i] = true
        #     mc.conf_down[site] = false
        #     ratio = getRatio(mc.Ham.U_down, U_downinvs, oldconfdown, site, i)
    end
    if ratio^2 < 1.0 && rand(ctx.rng) > ratio^2
        # measure acceptance rate here
        measure!(ctx, :acc, 0.0)
    else
        if flag == 1
            mc.conf_up[i] = false
            mc.conf_up[site] = true
            mc.conf_down[i] = true
            mc.conf_down[site] = false
        elseif flag == 2
            mc.conf_up[i] = true
            mc.conf_up[site] = false
            mc.conf_down[i] = false
            mc.conf_down[site] = true
        end
        measure!(ctx, :acc, 1.0)
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

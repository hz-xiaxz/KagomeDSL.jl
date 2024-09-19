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
    conf_up = FFS(rng, Ham.U_up)
    conf_down = FFS(rng, Ham.U_down)
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
    PBC = params[:PBC]
    lat = DoubleKagome(1.0, n1, n2, PBC)
    N_up = params[:N_up]
    N_down = params[:N_down]
    χ = params[:χ]
    Ham = Hamiltonian(χ, N_up, N_down, lat)
    ctx.rng = Random.Xoshiro(42)
    mc.conf_up = FFS(ctx.rng, Ham.U_up)
    mc.conf_down = FFS(ctx.rng, Ham.U_down)
end

@inline function Carlo.sweep!(mc::MC, ctx::MCContext)
    mc.conf_up = FFS(ctx.rng, mc.Ham.U_up)
    mc.conf_down = FFS(ctx.rng, mc.Ham.U_down)
    return nothing
end

@inline function Carlo.measure!(mc::MC, ctx::MCContext)
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

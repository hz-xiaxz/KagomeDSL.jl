mutable struct MC <: AbstractMC
    Ham::Hamiltonian
    conf::BitVector
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
    Ham = Hamiltonian(lat)
    rng = Random.Xoshiro(43)
    conf = FFS(rng, Ham.U)
    return MC(Ham, conf)
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
    PBC = params[:PBC]
    lat = DoubleKagome(1.0, n1, n2, PBC)
    Ham = Hamiltonian(lat)
    rng = Random.Xoshiro(43)
    ctx.rng = rng
    mc.conf = FFS(rng, Ham.U)
end

@inline function Carlo.sweep!(mc::MC, ctx::MCContext)
    mc.conf = FFS(ctx.rng, mc.Ham.U)
    return nothing
end

@inline function Carlo.measure!(mc::MC, ctx::MCContext)
    OL = getOL(mc.Ham, mc.conf)
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
    out["conf"] = Vector{Bool}(mc.conf)
    return nothing
end

@inline function Carlo.read_checkpoint!(mc::MC, in::HDF5.Group)
    mc.conf = BitVector(read(in, "conf"))
    return nothing
end

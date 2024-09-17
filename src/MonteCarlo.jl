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
    N_up = params[:N_up]
    N_down = params[:N_down]
    χ = params[:χ]
    Ham = Hamiltonian(χ, N_up, N_down, lat)
    rng = Random.Xoshiro(42)
    conf = vcat(FFS(rng, Ham.U_up), FFS(rng, Ham.U_down))
    return MC(Hamiltonian, conf)
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
    rng = Random.Xoshiro(42)
    ctx.rng = rng
    conf = vcat(FFS(rng, Ham.U_up), FFS(rng, Ham.U_down))
    mc.conf = conf
end

@inline function Carlo.sweep!(mc::MC{B}, ctx::MCContext) where {B}
    mc.conf = vcat(FFS(ctx.rng, mc.Ham.U_up), FFS(ctx.rng, mc.Ham.U_down))
    return nothing
end

@inline function Carlo.measure!(mc::MC{B}, ctx::MCContext) where {B}
    conf_up = FFS(ctx.rng, mc.Ham.U_up)
    conf_down = FFS(ctx.rng, mc.Ham.U_down)
    OL = getOL(mc.ham, conf_up, conf_down)
    measure!(ctx, :OL, OL)

    return nothing
end


@inline function Carlo.write_checkpoint(mc::MC{B}, out::HDF5.Group) where {B}
    out["conf"] = Vector{Bool}(mc.conf)
    return nothing
end

@inline function Carlo.read_checkpoint!(mc::MC{B}, in::HDF5.Group) where {B}
    mc.conf = BitVector(read(in, "conf"))
    return nothing
end

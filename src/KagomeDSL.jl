module KagomeDSL
using LinearAlgebra
using StatsBase
using Random
using Carlo
using HDF5
import LinearAlgebra.BLAS: geru!

export DoubleKagome
include("Lattice.jl")

export Hamiltonian, fast_update, Sz, spinInteraction!
include("Hamiltonian.jl")

export MC, MCContext, tilde_U
include("MonteCarlo.jl")
end

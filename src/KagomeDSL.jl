module KagomeDSL
using BitBasis
using LinearAlgebra
using StatsBase
using Random
using Carlo
using HDF5
export FFS
include("FFS.jl")

export Kagome, DoubleKagome, nearestNeighbor, AbstractLattice, ns
include("Lattice.jl")

export orbitals, Hamiltonian
include("Hamiltonian.jl")

export MC, MCContext
include("MonteCarlo.jl")
end

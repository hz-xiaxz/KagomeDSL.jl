module KagomeDSL
using MKL
using LinearAlgebra
using StatsBase
using Random
using Carlo
using HDF5
@info "Active BLAS vendor: $(BLAS.vendor())"
@info "BLAS.get_config(): $(BLAS.get_config())"
MKL.set_num_threads(1)

export DoubleKagome
include("Lattice.jl")

export Hamiltonian, fast_update, Sz, spinInteraction!
include("Hamiltonian.jl")

export MC, MCContext, tilde_U
include("MonteCarlo.jl")

export get_sublattice_indices, get_reciprocal_vectors, get_K_points, spin_structure_factor
include("Measurements.jl")
end

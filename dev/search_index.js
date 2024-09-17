var documenterSearchIndex = {"docs":
[{"location":"95-reference/#reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"95-reference/#Contents","page":"Reference","title":"Contents","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Modules = [KagomeDSL]","category":"page"},{"location":"95-reference/#KagomeDSL.DoubleKagome-Tuple{Float64, Int64, Int64, Tuple{Bool, Bool}}","page":"Reference","title":"KagomeDSL.DoubleKagome","text":"double triangle unit cell Kagome lattice\n\nnote n1 here is still the number of repititions of triangle in the a1 direction, so n1 is asserted to be even. The total number is n1 * n2 * 3 sites.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.MC-Tuple{AbstractDict}","page":"Reference","title":"KagomeDSL.MC","text":"MC(params::AbstractDict)\n\n\n\nCreate a Monte Carlo object\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#Carlo.init!-Tuple{MC, MCContext, AbstractDict}","page":"Reference","title":"Carlo.init!","text":"Carlo.init!(mc::MC, ctx::MCContext, params::AbstractDict)\n\n\n\nInitialize the Monte Carlo object params\n\nn1 : Int number of cells in x direction\nn2 : Int number of cells in y direction\nPBC : Tuple{Bool,2} boundary condition, e.g. (false, false)\nχ : Float64 hopping strength\nN_up : Int number of up spinons\nN_down : Int number of down spinons\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.FFS-Tuple{Random.AbstractRNG, AbstractMatrix}","page":"Reference","title":"KagomeDSL.FFS","text":"FFS([rng=default_rng()], U::AbstractMatrix)\n\nEmploying Fast Fermion Sampling Algorithm to sample free Fermions\n\nU : the sampling ensemble, a matrix of size L x N, where L is the number of energy states and N is the number of Fermions\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.Hmat-Tuple{DoubleKagome, Float64}","page":"Reference","title":"KagomeDSL.Hmat","text":"Hmat(lat::AbstractLattice, χ::Float64, N_up::Int, N_down::Int)\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.fast_update-Union{Tuple{T1}, Tuple{T}, Tuple{N}, Tuple{N1}, Tuple{D}, Tuple{AbstractMatrix, AbstractMatrix, Union{BitBasis.DitStr{D, N1, T1}, BitBasis.SubDitStr{D, N1, T1}}, BitBasis.BitStr{N, T}}} where {D, N1, N, T, T1}","page":"Reference","title":"KagomeDSL.fast_update","text":"fast_update(U::AbstractMatrix, Uinvs::AbstractMatrix, newconf::BitStr{N,T}, oldconf::BitStr{N,T}) where {N,T}\n\nFast computing technique from Becca and Sorella 2017\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.getOL-Tuple{Hamiltonian, BitVector, BitVector}","page":"Reference","title":"KagomeDSL.getOL","text":"The observable O_L = fracxHpsi_Gxpsi_G\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.getxprime-Union{Tuple{T}, Tuple{N}, Tuple{Hamiltonian, BitBasis.BitStr{N, T}}} where {N, T}","page":"Reference","title":"KagomeDSL.getxprime","text":"return x = Hx  where H = -t _ij χ_ij f_i^ f_j\n\n\n\n\n\n","category":"method"},{"location":"","page":"KagomeDSL","title":"KagomeDSL","text":"CurrentModule = KagomeDSL","category":"page"},{"location":"#KagomeDSL","page":"KagomeDSL","title":"KagomeDSL","text":"","category":"section"},{"location":"","page":"KagomeDSL","title":"KagomeDSL","text":"Documentation for KagomeDSL.","category":"page"}]
}

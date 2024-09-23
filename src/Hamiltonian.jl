# fluxed transition rule defined in the paper
# (i, j) i is the initial state, j is the final state
# in-cell transitions, needs index1 mod 6 == index2 mod 6
function bondNum(s1::Int, s2::Int)
    first_cell_num = (s1 - 1) ÷ 6
    label_s1 = (s1 - 1) % 6
    # if label_s1 <0 add 6 until it is positive
    while label_s1 < 0
        label_s1 += 6
    end
    new_s1 = label_s1 + 1
    new_s2 = s2 - 6 * first_cell_num
    return (new_s1, new_s2)
end

function setTunnel!(
    Tunneling::AbstractMatrix,
    s1::Int,
    s2::Int,
    PBCs1::Int,
    PBCs2::Int,
    link::Dict,
)
    bond_num = bondNum(PBCs1, PBCs2)
    if haskey(link, bond_num)
        # a bit hacky `=`
        Tunneling[s1, s2] = link[bond_num] * 1.0
    end

    bond_num2 = bondNum(PBCs2, PBCs1)
    if haskey(link, bond_num2)
        Tunneling[s2, s1] = 1.0 * link[bond_num2]
    end
end

"""
    Hmat(lat::AbstractLattice, χ::Float64, N_up::Int, N_down::Int)
"""
function Hmat(lat::DoubleKagome)

    n1 = lat.n1
    n2 = lat.n2
    ns = n1 * n2 * 3

    tunneling = zeros(Float64, ns, ns)

    link_in = Dict(
        (1, 2) => 1,
        (1, 3) => 1,
        (2, 3) => 1,
        (2, 4) => -1,
        (4, 6) => 1,
        (4, 5) => 1,
        (5, 6) => 1,
        (2, 1) => 1,
        (3, 1) => 1,
        (3, 2) => 1,
        (4, 2) => -1,
        (6, 4) => 1,
        (5, 4) => 1,
        (6, 5) => 1,
    )

    # maybe more elegent to modify bondNum function
    link_inter = Dict(
        (3, 3n1 + 1) => -1,
        (3, 3n1 - 1) => -1,
        (6, 3n1 + 2) => -1,
        (6, 3n1 + 4) => 1,
        (5, 7) => 1,
        (1, -1) => 1,
        (1, 3 - 3n1) => -1,
        (2, 6 - 3n1) => -1,
        (4, 6 - 3n1) => 1,
        (5, 9 - 3n1) => -1,
    )

    for bond in lat.nn
        s1, s2 = Tuple(bond)
        # consider in-cell cases
        if (s1 - 1) ÷ 6 == (s2 - 1) ÷ 6
            setTunnel!(tunneling, s1, s2, s1, s2, link_in)
        end
        # consider inter-cell cases
        # consider PBC

        # start from OBC

        if (s1 - 1) ÷ 6 != (s2 - 1) ÷ 6
            if lat.PBC == (false, false)
                setTunnel!(tunneling, s1, s2, s1, s2, link_inter)

            elseif lat.PBC == (true, false)
                for PBCs1 in (s1 - 3n1, s1 + 3n1, s1)
                    for PBCs2 in (s2 - 3n1, s2 + 3n1, s2)
                        setTunnel!(tunneling, s1, s2, PBCs1, PBCs2, link_inter)
                    end
                end

            elseif lat.PBC == (false, true)
                for PBCs1 in (s1 - 3ns, s1 + 3ns, s1)
                    for PBCs2 in (s2 - 3ns, s2 + 3ns, s2)
                        setTunnel!(tunneling, s1, s2, PBCs1, PBCs2, link_inter)
                    end
                end

            elseif lat.PBC == (true, true)
                for PBCs1 in (s1 - 3n1, s1 + 3n1, s1 - 3ns, s1 + 3ns, s1)
                    for PBCs2 in (s1 - 3n1, s1 + 3n1, s2 - 3ns, s2 + 3ns, s2)
                        setTunnel!(tunneling, s1, s2, PBCs1, PBCs2, link_inter)
                    end
                end
            end
        end
    end
    return tunneling
end

# temporarily separate the N_up and N_down subspaces
function orbitals(H_mat::Matrix{Float64})
    # get sampling ensemble U_up and U_down
    vals, vecs = eigen(H_mat)
    return vecs
end

# what's in the Hamiltonian
struct Hamiltonian
    U::Matrix{Float64}
    H_mat::Matrix{Float64}
    nn::AbstractArray
end

function Hamiltonian(lat::T) where {T<:AbstractLattice}
    H_mat = Hmat(lat)
    U = orbitals(H_mat)
    nn = lat.nn
    return Hamiltonian(U, H_mat, nn)
end


"""

return ``|x'> = H|x>``  where ``H = -t ∑_{<i,j>} χ_{ij} f_i^† f_j``
"""
@inline function getxprime(Ham::Hamiltonian, x::BitStr{N,T}) where {N,T}
    H_mat = Ham.H_mat
    nn = Ham.nn
    @assert N == 2 * size(H_mat)[1] "x should have the same 2x length as $(size(H_mat)[1]),  got: $N"
    L = length(x) ÷ 2  # Int division
    xprime = Dict{typeof(x),Float64}()
    # consider the spin up case
    @inbounds for i = 1:L÷2
        neigh_bond = filter(x -> x[1] == i, nn)
        for bond in neigh_bond
            if readbit(x, i) == 1 && readbit(x, bond[2]) == 0
                _x = x
                _x &= ~indicator(T, i)
                _x |= indicator(T, bond[2])
                xprime[_x] = get!(xprime, _x, 0.0) + H_mat[bond] # Hopping
            elseif readbit(x, i) == 0 && readbit(x, bond[2]) == 1
                _x = x
                _x &= ~indicator(T, bond[2])
                _x |= indicator(T, i)
                xprime[_x] = get!(xprime, _x, 0.0) + H_mat[bond[2], bond[1]] # Hopping
            end

            if readbit(x, i + L) == 1 && readbit(x, bond[2] + L) == 0
                _x = x
                _x &= ~indicator(T, i + L)
                _x |= indicator(T, bond[2] + L)
                xprime[_x] = get!(xprime, _x, 0.0) + H_mat[bond] # Hopping
            elseif readbit(x, i + L) == 0 && readbit(x, bond[2] + L) == 1
                _x = x
                _x &= ~indicator(T, bond[2] + L)
                _x |= indicator(T, i + L)
                xprime[_x] = get!(xprime, _x, 0.0) + H_mat[bond[2], bond[1]] # Hopping
            end
        end
    end
    return xprime
end

@inline function Gutzwiller(x::BitStr{N,T}) where {N,T}
    L = length(x) ÷ 2
    @inbounds for i = 1:L
        if readbit(x, i) == 1 && readbit(x, i + L) == 1
            return 0
        end
    end
    return 1
end

"""
    fast_update(U::AbstractMatrix, Uinvs::AbstractMatrix, newconf::BitStr{N,T}, oldconf::BitStr{N,T}) where {N,T}

Fast computing technique from Becca and Sorella 2017
"""
@inline function fast_update(
    U::AbstractMatrix,
    Uinvs::AbstractMatrix,
    newconf::Union{SubDitStr{D,N1,T1},DitStr{D,N1,T1}},
    oldconf::BitStr{N,T},
) where {D,N1,N,T,T1}
    @assert length(newconf) == N "The length of the new configuration should be the same as the old configuration, got: $(length(newconf))(old) and $N(new)"
    Rl = -1 # if not found should return error
    k = -1
    flag = 0
    @inbounds for i = 1:N
        if getindex(oldconf, i) == 1 && getindex(newconf, i) == 0
            Rl = i # the old position of the l-th electron
            flag += 1
        end
        if getindex(newconf, i) == 1 && getindex(oldconf, i) == 0
            k = i # the new position of the l-th electron, K = R_l'
            flag += 1
        end
        if flag == 2
            break
        end
    end
    if flag == 0
        return 1.0
    end
    l = sum(oldconf[1:Rl]) # l-th electron
    @show oldconf
    @show newconf
    ratio = sum(U[k, :] .* Uinvs[:, l])
    return ratio
end

@doc raw"""

The observable ``O_L = \frac{<x|H|\psi_G>}{<x|\psi_G>}``
"""
@inline function getOL(Ham::Hamiltonian, conf::BitVector)
    L = length(conf) ÷ 2
    @inbounds for i = 1:L
        if conf[i] == 1 && conf[i+L] == 1
            return 0.0
        end
    end
    conf = LongBitStr(vcat(conf_up, conf_down))
    OL = 0.0
    U_invs = Ham.U[conf, :] \ I # do invs more efficiently
    xprime = getxprime(Ham, conf)
    old_conf = LongBitStr(conf)
    @inbounds for (confstr, coff) in pairs(xprime)
        if Gutzwiller(confstr) == 0
            continue
        else
            OL += coff * fast_update(Ham.U, U_invs, confstr, conf)
        end
    end
    return OL
end

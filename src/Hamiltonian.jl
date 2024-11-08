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
    return -tunneling
    # Seems like the sign of tunneling is opposite to the one in the paper
end

# temporarily separate the N_up and N_down subspaces
function orbitals(H_mat::Matrix{Float64}, N_up::Int, N_down::Int)
    search_num = max(N_up, N_down)
    # get sampling ensemble U_up and U_down
    decomp, history =
        ArnoldiMethod.partialschur(H_mat, nev = search_num, tol = 1e-14, which = :SR)
    # select N lowest eigenvectors as the sampling ensemble
    U_up = decomp.Q[:, 1:N_up]
    U_down = decomp.Q[:, 1:N_down]
    return U_up, U_down
end

struct Hamiltonian
    N_up::Int
    N_down::Int
    U_up::Matrix{Float64}
    U_down::Matrix{Float64}
    H_mat::Matrix{Float64}
    nn::AbstractArray
end

function Hamiltonian(N_up::Int, N_down::Int, lat::T) where {T<:AbstractLattice}
    H_mat = Hmat(lat)
    U_up, U_down = orbitals(H_mat, N_up, N_down)
    nn = lat.nn
    return Hamiltonian(N_up, N_down, U_up, U_down, H_mat, nn)
end


# """
#     Splus(i::Int, x::BitStr{N,T}) where {N,T}
# ``S+ = f^†_{↑} f_{↓}``
# """
# function Splus(i::Int, x::BitStr{N,T}) where {N,T}
#     L = length(x) ÷ 2
#     if readbit(x, i) == 0 && readbit(x, i + L) == 1
#         _x = x
#         _x &= ~indicator(T, i + L)
#         _x |= indicator(T, i)
#         return _x, 1.0
#     else
#         return x, 0.0
#     end
# end

# """
#     Sminus(i::Int, x::BitStr{N,T}) where {N,T}
# ``S- = f^†_{↓} f_{↑}``
# """
# function Sminus(i::Int, x::BitStr{N,T}) where {N,T}
#     L = length(x) ÷ 2
#     if readbit(x, i) == 1 && readbit(x, i + L) == 0
#         return _x, 1.0
#     else
#         return x, 0.0
#     end
# end

"""
    Sz(i::Int, x::BitStr{N,T}) where {N,T}
``Sz = 1/2 (f^†_{↑} f_{↑} - f^†_{↓} f_{↓})``
"""
function Sz(i::Int, x::BitStr{N,T}) where {N,T}
    L = length(x) ÷ 2
    if readbit(x, i) == 1 && readbit(x, i + L) == 0
        return 1.0 / 2
    elseif readbit(x, i) == 0 && readbit(x, i + L) == 1
        return -1.0 / 2
    elseif readbit(x, i) == 1 && readbit(x, i + L) == 1
        error("site $i is doubly occupied")
        return 0.0
    else
        error("site $i is not occupied")
        return 0.0
    end
end

function SzInteraction!(xprime::Dict, x::BitStr{N,T}, i::Int, j::Int) where {N,T}
    xprime[(-1, -1, -1, -1)] = get!(xprime, (-1, -1, -1, -1), 0.0) + Sz(i, x) * Sz(j, x)
    return nothing
end

function spinInteraction!(xprime::Dict, x::BitStr{N,T}, i::Int, j::Int) where {N,T}
    L = length(x) ÷ 2
    # 1/2 (S+_i S-_j + S-_i S+_j)
    # first term
    # spin up at j, spin down at i
    if readbit(x, j) == 1 && readbit(x, i + L) == 1
        # K is the new index
        K_up = i
        K_down = j
        l_up = sum(x[1:j])
        l_down = sum(x[L+1:i+L])
        new_conf = (K_up, l_up, K_down, l_down)
        xprime[new_conf] = get!(xprime, new_conf, 0.0) + 1.0 / 2.0
    elseif readbit(x, j + L) == 1 && readbit(x, i) == 1
        K_up = j
        K_down = i
        l_up = sum(x[1:i])
        l_down = sum(x[L+1:j+L])
        new_conf = (K_up, l_up, K_down, l_down)
        xprime[new_conf] = get!(xprime, new_conf, 0.0) + 1.0 / 2.0
    end
    return nothing
end

"""

return ``|x'> = H|x>``  where ``H`` is the Heisenberg Hamiltonian
Note ``|x>`` here should also be a Mott state.
"""
@inline function getxprime(Ham::Hamiltonian, x::BitStr{N,T}) where {N,T}
    H_mat = Ham.H_mat
    nn = Ham.nn
    @assert N == 2 * size(H_mat)[1] "x should have the same 2x length as $(size(H_mat)[1]),  got: $N"
    xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()
    # just scan through all the bonds
    @inbounds for bond in nn
        @assert bond[2] > bond[1] "The second site should be larger than the first site, got: $(bond[2]) and $(bond[1])"
        spinInteraction!(xprime, x, bond[1], bond[2])
        spinInteraction!(xprime, x, bond[2], bond[1])
        SzInteraction!(xprime, x, bond[1], bond[2])
        SzInteraction!(xprime, x, bond[2], bond[1])
    end
    return xprime
end

@doc raw"""

The observable ``O_L = \frac{<x|H|\psi_G>}{<x|\psi_G>}``
The Hamiltonian should be the real one!
"""
@inline function getOL(mc::AbstractMC, conf_up::BitVector, conf_down::BitVector)
    @assert length(conf_up) == length(conf_down) "The length of the up and down configurations should be the same, got: $(length(conf_up)) and $(length(conf_down))"
    L = length(conf_up)
    # if double occupied state, no possibility to have a non-zero overlap
    # don't need to check this because MC will not propose double occupied states
    # any(conf_up .& conf_down) && return 0.0

    OL = 0.0
    oldconfstr = LongBitStr(vcat(conf_up, conf_down))
    xprime = getxprime(mc.Ham, oldconfstr)
    @inbounds for (conf, coff) in pairs(xprime)
        if conf == (-1, -1, -1, -1)
            OL += coff
        else
            # Gutzwiller(conf) == 0.0 && continue
            update_up = mc.W_up[conf[1], conf[2]]
            update_down = mc.W_down[conf[3], conf[4]]
            OL += coff * update_up * update_down
        end
    end
    return OL
end

# fluxed transition rule defined in the paper
# (i, j) i is the initial state, j is the final state
# in-cell transitions, needs index1 mod 6 == index2 mod 6
function bondNum(s1::Int, s2::Int)
    first_cell_num = (s1 - 1) ÷ 6
    if s1 <= 0 && mod(s1 - 1, 6) != 0
        first_cell_num -= 1
    end

    new_s1 = mod(s1 - 1, 6) + 1
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
    sign::Float64 = 1.0,
)
    bond_num = bondNum(PBCs1, PBCs2)
    if haskey(link, bond_num)
        Tunneling[s1, s2] = link[bond_num] * sign
    end

    bond_num2 = bondNum(PBCs2, PBCs1)
    if haskey(link, bond_num2)
        Tunneling[s2, s1] = link[bond_num2] * sign
    end
end

function get_boundary_shifts(lat::AbstractLattice, s1::Int, s2::Int, n1::Int, ns::Int)
    PBC1, PBC2 = lat.PBC
    anti1, anti2 = lat.antiPBC

    # No shifts for open boundary conditions
    if lat.PBC == (false, false)
        return [(s1, s2, 1.0)]
    end

    shifts = Tuple{Int,Int,Float64}[]

    # Handle PBC in first direction
    if PBC1
        base_shifts = [0, -3n1, 3n1]
        for shift in base_shifts
            sign = (anti1 && (shift != 0)) ? -1.0 : 1.0
            push!(shifts, (s1 + shift, s2, sign))
        end
    end

    # Handle PBC in second direction
    if PBC2
        base_shifts = [0, -3ns, 3ns]
        for shift in base_shifts
            sign = (anti2 && (shift != 0)) ? -1.0 : 1.0
            push!(shifts, (s1, s2 + shift, sign))
        end
    end

    # Handle both directions if both are periodic
    if PBC1 && PBC2
        for shift1 in [-3n1, 0, 3n1]
            for shift2 in [-3ns, 0, 3ns]
                sign = 1.0
                if anti1 && (shift1 != 0)
                    sign *= -1.0
                end
                if anti2 && (shift2 != 0)
                    sign *= -1.0
                end
                push!(shifts, (s1 + shift1, s2 + shift2, sign))
            end
        end
    end

    return unique(shifts)
end

function apply_boundary_conditions!(
    tunneling::AbstractMatrix,
    lat::AbstractLattice,
    s1::Int,
    s2::Int,
    n1::Int,
    ns::Int,
    link_in::Dict,
    link_inter::Dict,
)
    # Handle in-cell cases
    if (s1 - 1) ÷ 6 == (s2 - 1) ÷ 6
        setTunnel!(tunneling, s1, s2, s1, s2, link_in)
        return
    end

    # Handle inter-cell cases
    shifts = get_boundary_shifts(lat, s1, s2, n1, ns)
    # first check if s1 and s2 are linked, then do the boundary condition
    for (PBCs1, PBCs2, sign) in shifts
        if haskey(link_inter, bondNum(PBCs1, PBCs2))
            setTunnel!(tunneling, s1, s2, PBCs1, PBCs2, link_inter, sign)
        end
    end
end

"""
    Hmat(lat::DoubleKagome) -> Matrix{Float64}

Return the Spinor Hamiltonian matrix for a DoubleKagome lattice.
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
        apply_boundary_conditions!(tunneling, lat, s1, s2, n1, ns, link_in, link_inter)
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
    Sz(i::Int, kappa_up::Vector{Int}, kappa_down::Vector{Int}) -> Float64

``Sz = 1/2 (f^†_{↑} f_{↑} - f^†_{↓} f_{↓})``
Calculate the z-component of spin at site `i` given up and down spin configurations.

Returns:
    +1/2 for up spin
    -1/2 for down spin

Throws:
    ArgumentError: if site is doubly occupied or empty
    BoundsError: if i is outside the valid range
"""
@inline function Sz(i::Int, kappa_up::Vector{Int}, kappa_down::Vector{Int})
    # Bounds check
    n = length(kappa_up)
    @boundscheck begin
        1 ≤ i ≤ n || throw(BoundsError(kappa_up, i))
        length(kappa_down) == n || throw(
            DimensionMismatch(
                "kappa_up and kappa_down must have same length, got $n and $(length(kappa_down))",
            ),
        )
    end

    # Check occupation state
    up_occupied = is_occupied(kappa_up, i)
    down_occupied = is_occupied(kappa_down, i)

    # Use @fastmath for potential performance improvement in simple arithmetic
    return if up_occupied && !down_occupied
        0.5
    elseif !up_occupied && down_occupied
        -0.5
    elseif up_occupied && down_occupied
        throw(
            ArgumentError(
                "Site $i is doubly occupied, with kappa_up: $kappa_up and kappa_down: $kappa_down",
            ),
        )
    else
        throw(
            ArgumentError(
                "Site $i is unoccupied, with kappa_up: $kappa_up and kappa_down: $kappa_down",
            ),
        )
    end
end


function SzInteraction!(
    xprime::Dict,
    kappa_up::Vector{Int},
    kappa_down::Vector{Int},
    i::Int,
    j::Int,
)
    xprime[(-1, -1, -1, -1)] =
        get!(xprime, (-1, -1, -1, -1), 0.0) +
        Sz(i, kappa_up, kappa_down) * Sz(j, kappa_up, kappa_down)
    return nothing
end

"""
    spinInteraction!(xprime::Dict, kappa_up::Vector{Int}, kappa_down::Vector{Int}, i::Int, j::Int)

Compute the spin flip term 1/2(S+_i S-_j + S-_i S+_j) contribution to xprime.

Parameters:
- `xprime`: Dictionary to store the new configurations and their coefficients
- `kappa_up`: Configuration of up spins
- `kappa_down`: Configuration of down spins
- `i`, `j`: Sites where the spin interaction operates

The function handles two cases:
1. S+_i S-_j: when j has up spin and i has down spin
2. S-_i S+_j: when i has up spin and j has down spin

Each case contributes with coefficient 1/2.
"""
function spinInteraction!(
    xprime::Dict,
    kappa_up::Vector{Int},
    kappa_down::Vector{Int},
    i::Int,
    j::Int,
)
    # Case 1: S+_i S-_j
    # j has up spin (kappa_up[j] ≠ 0) and i has down spin (kappa_down[i] ≠ 0)
    if is_occupied(kappa_up, j) && is_occupied(kappa_down, i)
        # i, j are the original labels, in the {R_l} set
        # kappa[i], kappa[j] bookkeep the order inside tilde U, which is in the {l} set
        K_up = i    # i gets the up spin
        K_down = j  # j gets the down spin
        l_up = kappa_up[j]   # take the up index from j
        l_down = kappa_down[i]  # take the down index from i
        new_conf = (K_up, l_up, K_down, l_down)
        xprime[new_conf] = get!(xprime, new_conf, 0.0) + 1.0 / 2.0
    end

    # Case 2: S-_i S+_j
    # i has up spin (kappa_up[i] ≠ 0) and j has down spin (kappa_down[j] ≠ 0)
    if is_occupied(kappa_up, i) && is_occupied(kappa_down, j)
        K_up = j    # j gets the up spin
        K_down = i  # i gets the down spin
        l_up = kappa_up[i]   # take the up index from i
        l_down = kappa_down[j]  # take the down index from j
        new_conf = (K_up, l_up, K_down, l_down)
        xprime[new_conf] = get!(xprime, new_conf, 0.0) + 1.0 / 2.0
    end

    return nothing
end

"""

return ``|x'> = H|x>``  where ``H`` is the Heisenberg Hamiltonian
Note ``|x>`` here should also be a Mott state.
"""
@inline function getxprime(Ham::Hamiltonian, kappa_up::Vector{Int}, kappa_down::Vector{Int})
    nn = Ham.nn
    @assert length(kappa_up) == length(kappa_down) "The length of the up and down configurations should be the same, got: $(length(kappa_up)) and $(length(kappa_down))"
    xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()
    # just scan through all the bonds
    @inbounds for bond in nn
        @assert bond[2] > bond[1] "The second site should be larger than the first site, got: $(bond[2]) and $(bond[1])"
        spinInteraction!(xprime, kappa_up, kappa_down, bond[1], bond[2])
        spinInteraction!(xprime, kappa_up, kappa_down, bond[2], bond[1])
        SzInteraction!(xprime, kappa_up, kappa_down, bond[1], bond[2])
        SzInteraction!(xprime, kappa_up, kappa_down, bond[2], bond[1])
    end
    return xprime
end

@doc raw"""

The observable ``O_L = \frac{<x|H|\psi_G>}{<x|\psi_G>}``
The Hamiltonian should be the real one!
"""
@inline function getOL(mc::AbstractMC, kappa_up::Vector{Int}, kappa_down::Vector{Int})
    @assert length(kappa_up) == length(kappa_down) "The length of the up and down configurations should be the same, got: $(length(kappa_up)) and $(length(kappa_down))"
    # if double occupied state, no possibility to have a non-zero overlap
    # don't need to check this because MC will not propose double occupied states
    # any(conf_up .& conf_down) && return 0.0

    OL = 0.0
    xprime = getxprime(mc.Ham, kappa_up, kappa_down)
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

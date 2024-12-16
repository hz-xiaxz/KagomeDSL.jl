# fluxed transition rule defined in the paper
# (i, j) i is the initial state, j is the final state
# in-cell transitions, needs index1 mod 6 == index2 mod 6

function unitcell_coord(s::Int, n1::Int, n2::Int)
    @assert s >= 1 && s <= n1 * n2 * 6 "s should be in the range of 1 to ns, got: $s"
    unitcell_num = (s - 1) ÷ 6
    a1 = [4.0, 0.0]
    a2 = [1.0, sqrt(3.0)]
    unitcell_coord = (unitcell_num % n1) * a1 + (unitcell_num ÷ n1) * a2
    return unitcell_coord
end

function unitcell_diff(unitcell_coord1::Vector{Float64}, unitcell_coord2::Vector{Float64})
    diff = unitcell_coord1 - unitcell_coord2
    # projects onto a1 and a2
    a1 = [4.0, 0.0]
    a2 = [1.0, sqrt(3.0)]
    # Solve x*a1 + y*a2 = diff
    # [4.0 1.0  ] [x] = [diff[1]]
    # [0.0 √3.0 ] [y]   [diff[2]]

    # Solve using matrix inversion:
    # [x] = [4.0 1.0  ]^-1 [diff[1]]
    # [y]   [0.0 √3.0 ]    [diff[2]]

    dx = round(Int, (sqrt(3.0) * diff[1] - diff[2]) / (4.0 * sqrt(3.0)))
    dy = round(Int, diff[2] / sqrt(3.0))
    return dx, dy
end


function get_boundary_shifts(lat::AbstractLattice, s1::Int, s2::Int)

    PBC1, PBC2 = lat.PBC
    anti1, anti2 = lat.antiPBC
    # get unit cell coordinate

    n1 = lat.n1 ÷ 2
    n2 = lat.n2
    ns = n1 * n2 * 6
    @assert s1 >= 1 && s1 <= ns "s1 should be in the range of 1 to ns, got: $s1 in $ns"
    @assert s2 >= 1 && s2 <= ns "s2 should be in the range of 1 to ns, got: $s2 in $ns"
    u1 = unitcell_coord(s1, n1, n2)
    u2 = unitcell_coord(s2, n1, n2)
    dx, dy = unitcell_diff(u2, u1)

    # No shifts for open boundary conditions
    if lat.PBC == (false, false)
        return [(dx, dy, 1.0)]
    end

    shifts = [(dx, dy, 1.0)]
    # Handle PBC in first direction
    # Define shift ranges based on boundary conditions
    shifts_x = PBC1 ? [-n1, 0, n1] : [0]
    shifts_y = PBC2 ? [-n2, 0, n2] : [0]

    # Generate all combinations of shifts
    for shift1 in shifts_x, shift2 in shifts_y
        # Skip the no-shift case if it's already included
        (shift1 == 0 && shift2 == 0 && !isempty(shifts)) && continue

        # Calculate sign based on boundary crossings
        sign = 1.0
        if anti1 && shift1 != 0
            sign *= -1.0
        end
        if anti2 && shift2 != 0
            sign *= -1.0
        end

        push!(shifts, (dx + shift1, dy + shift2, sign))
    end
    # Check for inconsistent signs for same dx,dy pairs
    seen = Dict{Tuple{Int,Int},Float64}()
    for (dx, dy, sign) in shifts
        if haskey(seen, (dx, dy))
            if seen[(dx, dy)] != sign
                @warn "Inconsistent signs found for displacement ($dx,$dy)"
            end
        else
            seen[(dx, dy)] = sign
        end
    end
    return unique(shifts)
end

function apply_boundary_conditions!(
    tunneling::AbstractMatrix,
    lat::AbstractLattice,
    s1::Int,
    s2::Int,
    link_inter::Dict,
)
    cell1 = (s1 - 1) ÷ 6 + 1
    cell2 = (s2 - 1) ÷ 6 + 1
    # Handle in-cell cases
    @assert cell1 != cell2 "s1 and s2 should not be in the same cell, got: $cell1 and $cell2"
    label1 = (s1 - 1) % 6 + 1
    label2 = (s2 - 1) % 6 + 1
    shifts = get_boundary_shifts(lat, s1, s2)
    # check if the bond is linked
    for (dx, dy, sign) in shifts
        if haskey(link_inter, (label1, label2, dx, dy))
            tunneling[s1, s2] += sign * link_inter[(label1, label2, dx, dy)]
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
    # (lable1, label2, dx, dy)
    link_inter = Dict(
        (3, 5, -1, 1) => -1,
        (3, 1, 0, 1) => -1,
        (6, 2, 0, 1) => -1,
        (6, 4, 0, 1) => 1,
        (5, 1, 1, 0) => 1,
        (1, 5, -1, 0) => 1,
        (1, 3, 0, -1) => -1,
        (2, 6, 0, -1) => -1,
        (4, 6, 0, -1) => 1,
        (5, 3, 1, -1) => -1,
    )

    for cell1 = 1:(n1*n2÷2)
        sites1 = (cell1-1)*6+1:cell1*6
        for s1 in sites1, s2 in sites1
            s1 >= s2 && continue
            # s1 and s2 are in the same cell, so we only need to check link_in
            label1 = (s1 - 1) % 6 + 1
            label2 = (s2 - 1) % 6 + 1
            if haskey(link_in, (label1, label2))
                tunneling[s1, s2] = link_in[(label1, label2)]
            end
        end
    end
    for cell1 = 1:(n1*n2÷2)
        sites1 = (cell1-1)*6+1:cell1*6
        for cell2 = 1:(n1*n2÷2)
            sites2 = (cell2-1)*6+1:cell2*6
            cell1 == cell2 && continue
            for s1 in sites1, s2 in sites2
                s1 >= s2 && continue
                apply_boundary_conditions!(tunneling, lat, s1, s2, link_inter)
            end
        end
    end
    # careful!
    # verify tunneling matrix is upper triangular
    for i in axes(tunneling, 1)
        for j = 1:i-1
            if !iszero(tunneling[i, j])
                error("tunneling matrix must be upper triangular")
            end
        end
    end
    return -(tunneling + tunneling')
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

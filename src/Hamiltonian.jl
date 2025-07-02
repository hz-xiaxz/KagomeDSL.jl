# fluxed transition rule defined in the paper
# (i, j) i is the initial state, j is the final state
# in-cell transitions, needs index1 mod 6 == index2 mod 6

function unitcell_coord(lat::AbstractLattice, s::Int)
    n1 = lat.n1 ÷ 2
    n2 = lat.n2
    ns = n1 * n2 * 6
    @assert s >= 1 && s <= ns "s should be in the range of 1 to ns, got: $s"
    unitcell_num = (s - 1) ÷ 6
    a1 = lat.a1
    a2 = lat.a2
    unitcell_coord = (unitcell_num % n1) * a1 + (unitcell_num ÷ n1) * a2
    return unitcell_coord
end

function unitcell_diff(
    lat::AbstractLattice,
    unitcell_coord1::Vector{Float64},
    unitcell_coord2::Vector{Float64},
)
    diff = unitcell_coord1 - unitcell_coord2
    # projects onto a1 and a2
    a1 = lat.a1
    a2 = lat.a2
    # Solve x*a1 + y*a2 = diff
    # [a1[1] a2[1]] [x] = [diff[1]]
    # [a1[2] a2[2]] [y]   [diff[2]]

    # Solve using matrix inversion:
    # [x] = [a1[1] a2[1]]^-1 [diff[1]]
    # [y]   [a1[2] a2[2]]    [diff[2]]
    det = a1[1] * a2[2] - a1[2] * a2[1]
    dx = round(Int, (a2[2] * diff[1] - a2[1] * diff[2]) / det)
    dy = round(Int, (-a1[2] * diff[1] + a1[1] * diff[2]) / det)
    return dx, dy
end


function get_boundary_shifts(lat::AbstractLattice, s1::Int, s2::Int)
    @assert s1 != s2 "s1 and s2 should not be the same, got: $s1 and $s2"
    PBC1, PBC2 = lat.PBC
    anti1, anti2 = lat.antiPBC
    # get unit cell coordinate

    n1 = lat.n1 ÷ 2
    n2 = lat.n2
    ns = n1 * n2 * 6
    @assert s1 >= 1 && s1 <= ns "s1 should be in the range of 1 to ns, got: $s1 in $ns"
    @assert s2 >= 1 && s2 <= ns "s2 should be in the range of 1 to ns, got: $s2 in $ns"
    u1 = unitcell_coord(lat, s1)
    u2 = unitcell_coord(lat, s2)
    dx, dy = unitcell_diff(lat, u2, u1)

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
    B::Float64,
)
    cell1 = (s1 - 1) ÷ 6 + 1
    cell2 = (s2 - 1) ÷ 6 + 1
    @assert cell1 != cell2 "s1 and s2 should not be in the same cell"
    label1 = (s1 - 1) % 6 + 1
    label2 = (s2 - 1) % 6 + 1
    shifts = get_boundary_shifts(lat, s1, s2)

    r1 = get_site_coord(lat, s1)
    r_uc_1 = unitcell_coord(lat, s1)
    dr_2 = get_site_coord(lat, s2) - unitcell_coord(lat, s2)
    # note the phase is always calulated with one single bond
    for (dx, dy, sign) in shifts
        if haskey(link_inter, (label1, label2, dx, dy))
            r_uc_2_real = r_uc_1 + dx * lat.a1 + dy * lat.a2
            r2_real = r_uc_2_real + dr_2

            peierls_phase = (B / 2) * (r1[1] + r2_real[1]) * (r2_real[2] - r1[2])
            hopping_value = exp(im * peierls_phase)

            tunneling[s1, s2] += sign * link_inter[(label1, label2, dx, dy)] * hopping_value
        end
    end
end

const pi_link_in = Dict(
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
const pi_link_inter = Dict(
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

function get_site_coord(lat::AbstractLattice, s::Int)
    label = (s - 1) % 6
    uc_coord = unitcell_coord(lat, s)
    return uc_coord + lat.r[label+1]
end

"""
    Hmat(
        lat::DoubleKagome;
        link_in = pi_link_in,
        link_inter = pi_link_inter,
        B = 0.0,
    ) -> Matrix{ComplexF64}

Constructs the Spinon Hamiltonian matrix for a `DoubleKagome` lattice.

This function calculates the hopping terms within and between unit cells, incorporating
a Peierls phase to account for a magnetic field `B`. The resulting matrix
represents the Hamiltonian of the system.

# Arguments
- `lat::DoubleKagome`: The lattice structure for which to construct the Hamiltonian.
- `link_in`: A dictionary defining in-cell hopping terms. Defaults to `pi_link_in`.
- `link_inter`: A dictionary defining inter-cell hopping terms. Defaults to `pi_link_inter`.
- `B::Float64`: The magnetic field strength, used to calculate the Peierls phase. Defaults to `0.0`.

# Returns
- `Matrix{ComplexF64}`: The Hamiltonian matrix for the given lattice and parameters.
"""
function Hmat(lat::DoubleKagome; link_in = pi_link_in, link_inter = pi_link_inter, B = 0.0)
    n1 = lat.n1
    n2 = lat.n2
    ns = n1 * n2 * 3

    tunneling = zeros(ComplexF64, ns, ns)

    # in cell case
    for cell1 = 1:(n1*n2÷2)
        sites1 = ((cell1-1)*6+1):(cell1*6)
        for s1 in sites1, s2 in sites1
            s1 >= s2 && continue
            label1 = (s1 - 1) % 6 + 1
            label2 = (s2 - 1) % 6 + 1
            if haskey(link_in, (label1, label2))
                r1 = get_site_coord(lat, s1)
                r2 = get_site_coord(lat, s2)
                peierls_phase = (B / 2) * (r1[1] + r2[1]) * (r2[2] - r1[2])
                hopping_value = exp(im * peierls_phase)
                tunneling[s1, s2] = link_in[(label1, label2)] * hopping_value
            end
        end
    end

    # inter-cell case
    for cell1 = 1:(n1*n2÷2)
        sites1 = ((cell1-1)*6+1):(cell1*6)
        for cell2 = 1:(n1*n2÷2)
            sites2 = ((cell2-1)*6+1):(cell2*6)
            cell1 == cell2 && continue
            for s1 in sites1, s2 in sites2
                s1 >= s2 && continue
                apply_boundary_conditions!(tunneling, lat, s1, s2, link_inter, B)
            end
        end
    end
    # verify tunneling matrix is upper triangular
    for i in axes(tunneling, 1)
        for j = 1:(i-1)
            if !iszero(tunneling[i, j])
                error("tunneling matrix must be upper triangular")
            end
        end
    end
    @infiltrate
    return -(tunneling + tunneling')
    # from the sign of `-t`
end

# temporarily separate the N_up and N_down subspaces
function orbitals(H_mat::Matrix{ComplexF64}, N_up::Int, N_down::Int)
    search_num = max(N_up, N_down)
    # get sampling ensemble U_up and U_down
    decomp, history =
        ArnoldiMethod.partialschur(H_mat, nev = search_num, tol = 1e-14, which = :SR)
    # select N lowest eigenvectors as the sampling ensemble
    U_up = decomp.Q[:, 1:N_up]
    U_down = decomp.Q[:, 1:N_down]
    @infiltrate
    return U_up, U_down
end

struct Hamiltonian
    N_up::Int
    N_down::Int
    U_up::Matrix{ComplexF64}
    U_down::Matrix{ComplexF64}
    H_mat::Matrix{ComplexF64}
    nn::AbstractArray
end

function get_nn(H_mat::AbstractMatrix)
    # Get upper triangular non-zero elements
    indices = findall(!iszero, UpperTriangular(H_mat))
    return [(i[1], i[2]) for i in indices]
end


function Hamiltonian(
    N_up::Int,
    N_down::Int,
    lat::T;
    link_in = pi_link_in,
    link_inter = pi_link_inter,
    B = 0.0,
) where {T<:AbstractLattice}
    H_mat = Hmat(lat; link_in = link_in, link_inter = link_inter, B = B)
    U_up, U_down = orbitals(H_mat, N_up, N_down)
    nn = get_nn(H_mat)
    return Hamiltonian(N_up, N_down, U_up, U_down, H_mat, nn)
end


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
    xprime[ConfigKey(-1, -1, -1, -1)] =
        get!(xprime, ConfigKey(-1, -1, -1, -1), 0.0) +
        Sz(i, kappa_up, kappa_down) * Sz(j, kappa_up, kappa_down)
    return nothing
end

"""
    spinInteraction!(xprime::Dict, kappa_up::Vector{Int}, kappa_down::Vector{Int}, i::Int, j::Int)

Compute the spin flip term 1/2(S+_i S-_j + S-_i S+_j) contribution to xprime.


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
    i_up = @inbounds kappa_up[i]
    j_up = @inbounds kappa_up[j]
    i_down = @inbounds kappa_down[i]
    j_down = @inbounds kappa_down[j]
    # Case 1: S+_i S-_j
    # j has up spin (kappa_up[j] ≠ 0) and i has down spin (kappa_down[i] ≠ 0)
    if j_up != 0 && i_down != 0
        # i, j are the original labels, in the {R_l} set
        # kappa[i], kappa[j] bookkeep the order inside tilde U, which is in the {l} set
        # K_up = i    # i gets the up spin
        # K_down = j  # j gets the down spin
        # l_up = kappa_up[j]   # take the up index from j
        # l_down = kappa_down[i]  # take the down index from i
        new_conf = ConfigKey(i, j_up, j, i_down)
        xprime[new_conf] = get!(xprime, new_conf, 0.0) - 1.0 / 2.0
    end

    # Case 2: S-_i S+_j
    # i has up spin (kappa_up[i] ≠ 0) and j has down spin (kappa_down[j] ≠ 0)
    if i_up != 0 && j_down != 0
        # K_up = j    # j gets the up spin
        # K_down = i  # i gets the down spin
        # l_up = kappa_up[i]   # take the up index from i
        # l_down = kappa_down[j]  # take the down index from j
        new_conf = ConfigKey(j, i_up, i, j_down)
        xprime[new_conf] = get!(xprime, new_conf, 0.0) - 1.0 / 2.0
    end

    return nothing
end

struct ConfigKey
    K_up::Int
    l_up::Int
    K_down::Int
    l_down::Int
end
Base.hash(k::ConfigKey, h::UInt) =
    hash(k.l_down, hash(k.K_down, hash(k.l_up, hash(k.K_up, h))))
Base.:(==)(a::ConfigKey, b::ConfigKey) =
    a.K_up == b.K_up && a.l_up == b.l_up && a.K_down == b.K_down && a.l_down == b.l_down
Base.getindex(k::ConfigKey, i::Int) = getfield(k, i)
Base.length(::ConfigKey) = 4
Base.iterate(k::ConfigKey, state = 1) =
    state > 4 ? nothing : (getfield(k, state), state + 1)

"""

return ``|x'> = H|x>``  where ``H`` is the Heisenberg Hamiltonian
Note ``|x>`` here should also be a Mott state.
"""
@inline function getxprime(Ham::Hamiltonian, kappa_up::Vector{Int}, kappa_down::Vector{Int})
    nn = Ham.nn
    xprime = Dict{ConfigKey,Float64}()
    # just scan through all the bonds
    @inbounds for bond in nn
        spinInteraction!(xprime, kappa_up, kappa_down, bond[1], bond[2])
        SzInteraction!(xprime, kappa_up, kappa_down, bond[1], bond[2])
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

    OL = 0.0
    xprime = getxprime(mc.Ham, kappa_up, kappa_down)
    @inbounds for (conf, coff) in pairs(xprime)
        if conf == ConfigKey(-1, -1, -1, -1)
            OL += coff
        else
            update_up = mc.W_up[conf[1], conf[2]]
            update_down = mc.W_down[conf[3], conf[4]]
            OL += coff * update_up * update_down
        end
    end
    return OL
end

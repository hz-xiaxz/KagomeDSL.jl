# Spinon Hamiltonian for DoubleKagome lattice with π-flux per plaquette
# This implements the theoretical model from quantum spin liquid research
# where spinons (fractionalized spins) hop on the Kagome lattice with specific flux patterns.
# 
# Key concepts for AI understanding:
# - Spinons: Fractionalized spin excitations that emerge in quantum spin liquids
# - π-flux: Phase factor of π acquired when a spinon completes a loop around a plaquette
# - Unit cell: Contains 6 sites arranged in two triangular motifs
# - Peierls phase: Additional phase factor from magnetic field B

"""    
    unitcell_coord(lat::AbstractLattice, s::Int) -> Vector{Float64}

Compute the real-space coordinate of the unit cell containing site `s`.

For DoubleKagome lattice:
- Sites 1-6 belong to unit cell 1, sites 7-12 to unit cell 2, etc.
- Each unit cell has 6 sites arranged in two triangular motifs
- Returns the coordinate of the unit cell origin in real space

# Arguments
- `lat::AbstractLattice`: The lattice structure
- `s::Int`: Site index (1-indexed)

# Returns
- `Vector{Float64}`: Real-space coordinate of the unit cell origin
"""
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

"""    
    unitcell_diff(lat, unitcell_coord1, unitcell_coord2) -> (Int, Int)

Calculate the lattice vector difference between two unit cells in terms of basis vectors (a1, a2).

This function solves the linear system: unitcell_coord1 - unitcell_coord2 = dx*a1 + dy*a2
to find integer coefficients (dx, dy) representing the separation in lattice units.

# Arguments
- `lat::AbstractLattice`: The lattice structure containing basis vectors a1, a2
- `unitcell_coord1::Vector{Float64}`: Real-space coordinate of first unit cell
- `unitcell_coord2::Vector{Float64}`: Real-space coordinate of second unit cell

# Returns
- `(dx::Int, dy::Int)`: Lattice vector coefficients such that coord1 - coord2 = dx*a1 + dy*a2
"""
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


"""    
    get_boundary_shifts(lat::AbstractLattice, s1::Int, s2::Int) -> Vector{Tuple{Int,Int,Float64}}

Calculate all possible lattice vector shifts between sites s1 and s2 under periodic/antiperiodic boundary conditions.

For quantum spin systems with twisted boundary conditions:
- Periodic BC: ψ(r + L) = ψ(r)
- Antiperiodic BC: ψ(r + L) = -ψ(r)

Returns all equivalent separations considering lattice periodicity, with associated phase factors.

# Arguments
- `lat::AbstractLattice`: Lattice with boundary condition specifications
- `s1::Int, s2::Int`: Site indices for the hopping term

# Returns
- `Vector{Tuple{Int,Int,Float64}}`: List of (dx, dy, sign) tuples where:
  - `dx, dy`: Lattice vector coefficients
  - `sign`: Phase factor (±1) from antiperiodic boundary crossings
"""
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

"""    
    apply_boundary_conditions!(tunneling, lat, s1, s2, link_inter, B)

Apply boundary conditions and magnetic field effects to inter-cell hopping terms.

This function modifies the tunneling matrix by adding all possible hopping paths 
between sites s1 and s2, considering:
1. Periodic/antiperiodic boundary conditions
2. Peierls phase factors from magnetic field B
3. π-flux pattern encoded in link_inter dictionary

# Arguments
- `tunneling::AbstractMatrix`: Hopping matrix to be modified (in-place)
- `lat::AbstractLattice`: Lattice structure with boundary conditions
- `s1::Int, s2::Int`: Source and target site indices
- `link_inter::Dict`: Inter-cell hopping amplitudes with (label1, label2, dx, dy) keys
- `B::Float64`: Magnetic field strength for Peierls phase calculation

# Physics Details
- Peierls phase: exp(iB/2 * (x1+x2)(y2-y1)) accounts for vector potential A = (0, Bx)
- Each boundary crossing can contribute ±1 phase factor from antiperiodic BC
"""
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

# π-flux intra-cell hopping amplitudes for DoubleKagome lattice
# Each unit cell has 6 sites arranged as two triangular motifs:
# Sites 1,2,3 form one triangle; sites 4,5,6 form another
# The hopping pattern creates π-flux through specific plaquettes
#
# Key for AI understanding:
# - Positive values (+1): Normal hopping amplitude
# - Negative values (-1): Hopping with π phase (sign flip)
# - The pattern ensures each triangular plaquette accumulates π flux
# - Symmetric entries (i,j) and (j,i) have same amplitude (Hermitian hopping)
const pi_link_in = Dict(
    (1, 2) => 1,   # Triangle 1: sites 1-2-3
    (1, 3) => 1,
    (2, 3) => 1,
    (2, 4) => -1,  # Inter-triangle connection with π phase
    (4, 6) => 1,   # Triangle 2: sites 4-5-6  
    (4, 5) => 1,
    (5, 6) => 1,
    # Hermitian conjugates
    (2, 1) => 1,
    (3, 1) => 1,
    (3, 2) => 1,
    (4, 2) => -1,
    (6, 4) => 1,
    (5, 4) => 1,
    (6, 5) => 1,
)

# π-flux inter-cell hopping amplitudes for DoubleKagome lattice
# Format: (site_label1, site_label2, dx, dy) => hopping_amplitude
# where (dx, dy) are lattice vector coefficients: target_cell = source_cell + dx*a1 + dy*a2
#
# Physical interpretation:
# - These hoppings connect sites in different unit cells
# - The π-flux constraint requires specific sign patterns
# - Negative amplitudes (-1) create π phase shifts
# - The pattern ensures flux conservation across the entire lattice
#
# Example: (3, 5, -1, 1) => -1 means:
# - Hopping from site 3 to site 5 in the cell at (-a1 + a2) relative position
# - Amplitude is -1 (includes π phase)
const pi_link_inter = Dict(
    (3, 5, -1, 1) => -1,  # Site 3 to site 5 in neighboring cell
    (3, 1, 0, 1) => -1,   # Site 3 to site 1 in cell at +a2
    (6, 2, 0, 1) => -1,   # Site 6 to site 2 in cell at +a2
    (6, 4, 0, 1) => 1,    # Site 6 to site 4 in cell at +a2
    (5, 1, 1, 0) => 1,    # Site 5 to site 1 in cell at +a1
    (1, 5, -1, 0) => 1,   # Site 1 to site 5 in cell at -a1
    (1, 3, 0, -1) => -1,  # Site 1 to site 3 in cell at -a2
    (2, 6, 0, -1) => -1,  # Site 2 to site 6 in cell at -a2
    (4, 6, 0, -1) => 1,   # Site 4 to site 6 in cell at -a2
    (5, 3, 1, -1) => -1,  # Site 5 to site 3 in cell at +a1-a2
)

"""    
    get_site_coord(lat::AbstractLattice, s::Int) -> Vector{Float64}

Calculate the real-space coordinate of site `s` in the lattice.

Combines the unit cell position with the intra-cell site offset to give
the absolute position in real space. Essential for calculating Peierls phases
and analyzing spatial correlations.

# Arguments
- `lat::AbstractLattice`: Lattice structure
- `s::Int`: Site index (1-indexed)

# Returns
- `Vector{Float64}`: [x, y] coordinate in real space
"""
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
    return -(tunneling + tunneling')
    # from the sign of `-t`
end

"""    
    orbitals(H_mat::Matrix{ComplexF64}, N_up::Int, N_down::Int) -> (Matrix{ComplexF64}, Matrix{ComplexF64})

Compute the occupied spinon orbitals for the quantum spin liquid ground state.

In the spinon mean-field theory:
1. Spins are fractionalized into spinons (fermionic particles carrying spin-1/2)
2. The Hamiltonian H_mat describes spinon hopping on the lattice
3. Ground state is a filled Fermi sea of the lowest-energy spinon states
4. Separate up and down spinon sectors (SU(2) symmetry)

# Arguments
- `H_mat::Matrix{ComplexF64}`: Single-particle spinon Hamiltonian matrix
- `N_up::Int`: Number of up-spin spinons to fill
- `N_down::Int`: Number of down-spin spinons to fill

# Returns
- `(U_up, U_down)`: Matrices whose columns are the occupied spinon orbitals
  - `U_up`: N_up lowest eigenvectors for up spinons
  - `U_down`: N_down lowest eigenvectors for down spinons

# Physics Notes
- The spinon filling typically corresponds to the original spin-1/2 density
- For a spin-1/2 system: N_up + N_down = total number of spins
- The specific choice of N_up, N_down determines the magnetic properties
"""
function orbitals(H_mat::Matrix{ComplexF64}, N_up::Int, N_down::Int)
    # search_num = max(N_up, N_down)
    # get sampling ensemble U_up and U_down
    F = eigen(Hermitian(H_mat))
    p = sortperm(F.values)
    evalues = F.values[p]
    evecs = F.vectors[:, p]
    # select N lowest eigenvectors as the sampling ensemble
    U_up = evecs[:, 1:N_up]
    U_down = evecs[:, 1:N_down]
    return U_up, U_down
end

"""    
    Hamiltonian

Complete specification of the spinon mean-field Hamiltonian for quantum Monte Carlo simulations.

This structure contains all information needed for Variational 
Monte Carlo calculations of quantum spin liquid properties:

# Fields
- `N_up::Int`: Number of up-spin spinons (determines spin sector)
- `N_down::Int`: Number of down-spin spinons
- `U_up::Matrix{ComplexF64}`: Occupied up-spinon orbitals (columns are eigenvectors)
- `U_down::Matrix{ComplexF64}`: Occupied down-spinon orbitals
- `H_mat::Matrix{ComplexF64}`: Single-particle spinon Hamiltonian matrix
- `nn::AbstractArray`: List of nearest-neighbor bonds for efficient iteration

# Usage in Monte Carlo
- U_up, U_down define the reference state (filled Fermi sea)
- H_mat provides the hopping amplitudes for Monte Carlo updates
- nn specifies which bonds to consider in interaction terms

# Physical Interpretation
- Represents a quantum spin liquid state where spins are fractionalized
- Correlations arise from quantum fluctuations around the mean-field state
"""
struct Hamiltonian{N_up,N_down}
    U_up::Matrix{ComplexF64}
    U_down::Matrix{ComplexF64}
    H_mat::Matrix{ComplexF64}
    nn::AbstractArray
    U_up_plus::Matrix{ComplexF64}
    U_down_minus::Matrix{ComplexF64}
end

"""    
    get_nn(H_mat::AbstractMatrix) -> Vector{Tuple{Int,Int}}

Extract nearest-neighbor bond list from the Hamiltonian matrix.

Finds all non-zero off-diagonal elements in the upper triangular part of H_mat,
which correspond to hopping terms between connected sites. Essential for
efficient Monte Carlo sampling since we only need to consider active bonds.

# Arguments
- `H_mat::AbstractMatrix`: Hamiltonian matrix (typically hermitian)

# Returns
- `Vector{Tuple{Int,Int}}`: List of (i,j) pairs where i < j and H_mat[i,j] ≠ 0

# Notes
- Only upper triangular elements to avoid double-counting
"""
function get_nn(H_mat::AbstractMatrix)
    # Get upper triangular non-zero elements
    nn = Tuple{Int,Int}[]
    for j in axes(H_mat, 2)
        for i = 1:(j-1)
            if !iszero(H_mat[i, j])
                push!(nn, (i, j))
            end
        end
    end
    return nn
end

abstract type AbstractOperator end

struct SpinPlusOperator{N_up,N_down} <: AbstractOperator
    site::Int
end


"""    
    Hamiltonian(N_up, N_down, lat; link_in=pi_link_in, link_inter=pi_link_inter, B=0.0)

Construct a complete Hamiltonian structure for quantum Monte Carlo simulations.

This is the main constructor that builds everything needed for SSE Monte Carlo:
1. Constructs the single-particle spinon Hamiltonian matrix
2. Diagonalizes it to find the occupied orbitals
3. Extracts the nearest-neighbor bond structure
4. Packages everything for efficient Monte Carlo usage

# Arguments
- `N_up::Int, N_down::Int`: Number of up/down spinons (determines magnetic sector)
- `lat::AbstractLattice`: Lattice structure (typically DoubleKagome)
- `link_in`: Intra-cell hopping dictionary (default: π-flux pattern)
- `link_inter`: Inter-cell hopping dictionary (default: π-flux pattern)
- `B::Float64`: Magnetic field strength for Peierls phases (default: 0.0)

# Returns
- `Hamiltonian`: Complete structure ready for Monte Carlo simulations

# Physics Notes
- The choice of N_up, N_down determines the spin sector being studied
- For spin-1/2 systems: N_up + N_down = number of original spins
- Different (N_up, N_down) can access different quantum phases
"""
function Hamiltonian(
    N_up::Int,
    N_down::Int,
    lat::AbstractLattice;
    link_in = pi_link_in,
    link_inter = pi_link_inter,
    B = 0.0,
)
    return Hamiltonian(
        Val(N_up),
        Val(N_down),
        lat;
        link_in = link_in,
        link_inter = link_inter,
        B = B,
    )
end

function Hamiltonian(
    ::Val{N_up},
    ::Val{N_down},
    lat::AbstractLattice;
    link_in = pi_link_in,
    link_inter = pi_link_inter,
    B = 0.0,
) where {N_up,N_down}
    H_mat = Hmat(lat; link_in = link_in, link_inter = link_inter, B = B)
    U_up, U_down = orbitals(H_mat, N_up, N_down)

    ns = size(H_mat, 1)
    # Pre-compute orbitals for neighboring sectors needed for S+ transitions
    # U_up_plus: orbitals for N_up+1 sector (needed when S+ increases up-spin count)
    # U_down_minus: orbitals for N_down-1 sector (needed when S+ decreases down-spin count)
    if N_down > 0
        U_up_plus, U_down_minus = orbitals(H_mat, N_up + 1, N_down - 1)
    else
        # If N_down = 0, we can't have S+ transitions (no down spins to flip)
        U_up_plus = zeros(ComplexF64, ns, N_up + 1)
        U_down_minus = zeros(ComplexF64, ns, 0)  # N_down - 1 = -1, invalid
    end

    nn = get_nn(H_mat)
    return Hamiltonian{N_up,N_down}(U_up, U_down, H_mat, nn, U_up_plus, U_down_minus)
end


"""    
    Sz(i::Int, kappa_up::Vector{Int}, kappa_down::Vector{Int}) -> Float64

Calculate the z-component of spin at site `i` in the spinon representation.

In the spinon formulation, the original spin operators are written as:
- S^z_i = 1/2 (f^†_{i↑} f_{i↑} - f^†_{i↓} f_{i↓})
- S^+_i = f^†_{i↑} f_{i↓}, S^-_i = f^†_{i↓} f_{i↑}

where f_{iσ} are spinon annihilation operators. The constraint is exactly one 
spinon per site: n_{i↑} + n_{i↓} = 1 (no double occupancy or empty sites).

# Arguments
- `i::Int`: Site index (1-indexed)
- `kappa_up::Vector{Int}`: Up-spinon configuration (0 = empty, nonzero = occupied)
- `kappa_down::Vector{Int}`: Down-spinon configuration

# Returns
- `Float64`: +0.5 for up spin, -0.5 for down spin

# Exceptions
- `ArgumentError`: If site is doubly occupied or empty (violates spinon constraint)
- `BoundsError`: If i is outside valid range [1, length(kappa_up)]
- `DimensionMismatch`: If kappa_up and kappa_down have different lengths

# Physics Notes
- This enforces the constraint that each site has exactly one spinon
- Configurations violating the constraint (empty or doubly occupied) are unphysical
- Used in Monte Carlo to measure local magnetization and Ising interactions
- The spinon constraint is fundamental to the validity of the spin liquid state

# Performance Notes
- Marked @inline for performance in Monte Carlo loops
- Uses @boundscheck for optional bounds checking
- Calls is_occupied() helper function for spinon occupancy detection
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


"""    
    SzInteraction!(xprime, kappa_up, kappa_down, i, j)

Compute the Ising (S^z S^z) interaction term for the Heisenberg Hamiltonian.

For the Heisenberg model H = J ∑_{<i,j>} (S^+_i S^-_j + S^-_i S^+_j + S^z_i S^z_j),
this function handles the S^z_i S^z_j term, which is diagonal in the spinon basis.

The Ising term doesn't change the spinon configuration, so it contributes to the
identity component of the operator expansion (stored with key (-1,-1,-1,-1)).

# Arguments  
- `xprime::Dict`: Dictionary storing operator expansion coefficients
- `kappa_up::Vector{Int}`: Current up-spinon configuration
- `kappa_down::Vector{Int}`: Current down-spinon configuration  
- `i::Int, j::Int`: Sites for the interaction

# Side Effects
- Modifies xprime[(-1,-1,-1,-1)] by adding Sz(i) * Sz(j)

# Physics Notes
- This is the "easy-axis" part of the Heisenberg interaction
- Diagonal in the occupation basis, doesn't create/destroy spinons
- Combined with spinInteraction!() to give the full Heisenberg model
"""
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
    spinInteraction!(xprime, kappa_up, kappa_down, i, j)

Compute the transverse (spin-flip) part of the Heisenberg interaction.

For the Heisenberg Hamiltonian H = J ∑_{<i,j>} (S^+_i S^-_j + S^-_i S^+_j + S^z_i S^z_j),
this function handles the transverse terms S^+_i S^-_j + S^-_i S^+_j.

In the spinon representation:
- S^+_i S^-_j = f^†_{i↑} f_{i↓} f^†_{j↓} f_{j↑} (flips spins at both sites)
- S^-_i S^+_j = f^†_{i↓} f_{i↑} f^†_{j↑} f_{j↓} (flips spins at both sites)

These terms change the spinon configuration and are the source of quantum fluctuations.

# Arguments
- `xprime::Dict`: Dictionary storing new configurations and their amplitudes
- `kappa_up::Vector{Int}`: Current up-spinon configuration
- `kappa_down::Vector{Int}`: Current down-spinon configuration
- `i::Int, j::Int`: Sites for the spin-flip interaction

# Side Effects
- Adds entries to xprime for each allowed spin-flip process
- Key format: (new_up_site, old_up_orbital, new_down_site, old_down_orbital)
- Amplitude: -1/2 for each allowed process (negative from Heisenberg coupling)

# Physics Notes
- Only processes that respect the constraint (one spinon per site) are allowed
- Creates quantum entanglement between different spinon configurations
- Essential for accessing quantum spin liquid physics beyond mean-field
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
        new_conf = (i, j_up, j, i_down)
        xprime[new_conf] = get!(xprime, new_conf, 0.0) - 1.0 / 2.0
    end

    # Case 2: S-_i S+_j
    # i has up spin (kappa_up[i] ≠ 0) and j has down spin (kappa_down[j] ≠ 0)
    if i_up != 0 && j_down != 0
        # K_up = j    # j gets the up spin
        # K_down = i  # i gets the down spin
        # l_up = kappa_up[i]   # take the up index from i
        # l_down = kappa_down[j]  # take the down index from j
        new_conf = (j, i_up, i, j_down)
        xprime[new_conf] = get!(xprime, new_conf, 0.0) - 1.0 / 2.0
    end

    return nothing
end


"""    
    getxprime(Ham::Hamiltonian, kappa_up, kappa_down) -> Dict{NTuple{4,Int}, Float64}

Compute the action of the Heisenberg Hamiltonian on a spinon configuration.

This is the core function for Stochastic Series Expansion (SSE) Monte Carlo.
Given a spinon configuration |κ_up, κ_down⟩, it computes all configurations 
|κ'_up, κ'_down⟩ that can be reached by applying H, along with their amplitudes.

The result is H|κ⟩ = Σ_κ' xprime[κ'] |κ'⟩, where xprime encodes the expansion.

# Arguments
- `Ham::Hamiltonian`: Complete Hamiltonian specification
- `kappa_up::Vector{Int}`: Current up-spinon configuration
- `kappa_down::Vector{Int}`: Current down-spinon configuration

# Returns
- `Dict{NTuple{4,Int}, Float64}`: Dictionary mapping configurations to amplitudes
  - Key (-1,-1,-1,-1): Diagonal contribution (Ising terms)
  - Key (i,l_up,j,l_down): Off-diagonal contribution from spin-flip at sites i,j

# Algorithm
1. Iterate over all nearest-neighbor bonds in Ham.nn
2. For each bond, compute Ising (S^z S^z) and transverse (S^+ S^- + S^- S^+) terms
3. Accumulate results in xprime dictionary

# Usage in Monte Carlo
- Called during each Monte Carlo update to propose new configurations
- The amplitudes determine acceptance probabilities for updates
- Essential for sampling the quantum mechanical evolution

# Physics Notes
- |x⟩ must be a valid Mott state (one spinon per site)
- The Hamiltonian H is the original Heisenberg model, not the mean-field one
- Combines both diagonal (Ising) and off-diagonal (transverse) contributions
"""
@inline function getxprime(Ham::Hamiltonian, kappa_up::Vector{Int}, kappa_down::Vector{Int})
    nn = Ham.nn
    xprime = Dict{NTuple{4,Int},Float64}()
    # just scan through all the bonds
    @inbounds for bond in nn
        spinInteraction!(xprime, kappa_up, kappa_down, bond[1], bond[2])
        SzInteraction!(xprime, kappa_up, kappa_down, bond[1], bond[2])
    end
    return xprime
end


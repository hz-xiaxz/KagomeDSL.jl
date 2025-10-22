abstract type AbstractLattice end

function validate_boundary_conditions(PBC::Tuple{Bool,Bool}, antiPBC::Tuple{Bool,Bool})
    for (i, (pbc, apbc)) in enumerate(zip(PBC, antiPBC))
        if apbc && !pbc
            dir = i == 1 ? "first" : "second"
            throw(
                ArgumentError(
                    "Invalid boundary conditions in $dir direction: " *
                    "Cannot have antiperiodic boundary conditions without periodic boundary conditions",
                ),
            )
        end
    end
end

struct DoubleKagome <: AbstractLattice
    # number of repititions in directions of a1 and a2
    n1::Int # a1
    n2::Int # a2
    # parameter that defines the equilateral triangle side length
    t::Float64
    # translation vectors
    a1::Array{Float64,1}
    a2::Array{Float64,1}
    # coordinates of sites inside a unit cell
    r::Array{Array{Float64,1},1}

    PBC::Tuple{Bool,Bool}
    antiPBC::Tuple{Bool,Bool}
end

"""
    DoubleKagome(t::Float64, n1::Int, n2::Int, PBC::Tuple{Bool,Bool}; antiPBC=(false,false), trunc=Inf)

Construct a double unit cell Kagome lattice with 6 sites per unit cell.

Unlike a single unit cell Kagome lattice (which has 3 sites per unit cell forming one triangle),
the DoubleKagome lattice contains 6 sites per unit cell arranged in two triangular sublattices.
This structure is useful for studying systems with enlarged unit cells or specific magnetic
ordering patterns that require the doubled cell geometry.

# Arguments
- `t::Float64`: Parameter defining the equilateral triangle side length (lattice constant = 2t)
- `n1::Int`: Number of unit cell repetitions in the a1 direction (must be even)
- `n2::Int`: Number of unit cell repetitions in the a2 direction  
- `PBC::Tuple{Bool,Bool}`: Periodic boundary conditions in (a1, a2) directions
- `antiPBC::Tuple{Bool,Bool}`: Antiperiodic boundary conditions (default: (false, false))
- `trunc::Float64`: Truncation parameter (default: Inf, unused in current implementation)

# Returns
- `DoubleKagome`: Lattice structure with `n1 * n2 * 3` total sites

# Notes
- The constraint `n1 % 2 == 0` ensures proper tiling of the double unit cell
- Total number of sites is `n1 * n2 * 3` (though each unit cell has 6 sites, `n1` counts double-sized cells)
- Actual number of unit cells is `(n1 ÷ 2) * n2`, each containing 6 sites
- Lattice vectors: a1 = [4t, 0], a2 = [t, √3*t] 
- Each unit cell contains 6 sites arranged in two triangular motifs

# Example
```julia
# Create a 4×3 DoubleKagome lattice with periodic boundaries
lat = DoubleKagome(1.0, 4, 3, (true, true))
# Total sites: 4 * 3 * 3 = 36 sites (2 unit cells × 3 repetitions × 6 sites per unit cell)
```
"""
function DoubleKagome(
    t::Float64,
    n1::Int,
    n2::Int,
    PBC::Tuple{Bool,Bool};
    antiPBC::Tuple{Bool,Bool} = (false, false),
    trunc::Float64 = Inf,
)
    validate_boundary_conditions(PBC, antiPBC)
    @assert n1 % 2 == 0 "n1 must be even in DoubleKagome"
    a = 2t
    a1 = [2a, 0.0]
    a2 = [0.5a, 0.5*√3a]

    # coordinates of each site in the unit cell
    r1 = [0.0, 0.0]
    r2 = 0.25 * a1
    r3 = 0.5 * a2
    r4 = 0.5 * a1
    r5 = 0.75 * a1
    r6 = [2.5t, 0.5*√3t]
    r = [r1, r2, r3, r4, r5, r6]

    return DoubleKagome(n1, n2, t, a1, a2, r, PBC, antiPBC)
end


# For DoubleKagome: 6 sites per unit cell, but n1 counts double-sized cells
# so total sites = (n1 ÷ 2) * n2 * 6 = n1 * n2 * 3
ns(lat::DoubleKagome) = lat.n1 * lat.n2 * 3

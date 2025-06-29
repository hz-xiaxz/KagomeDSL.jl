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
double triangle unit cell Kagome lattice

note `n1` here is still the number of repititions of triangle in the `a1` direction, so `n1` is asserted to be even. The total number is `n1 * n2 * 3` sites.
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
    a2 = [0.5a, 0.5√3a]

    # coordinates of each site in the unit cell
    r1 = [0.0, 0.0]
    r2 = 0.25 * a1
    r3 = 0.5 * a2
    r4 = 0.5 * a1
    r5 = 0.75 * a1
    r6 = [2.5t, 0.5√3t]
    r = [r1, r2, r3, r4, r5, r6]

    return DoubleKagome(n1, n2, t, a1, a2, r, PBC, antiPBC)
end


ns(lat::AbstractLattice) = lat.n1 * lat.n2 * 3

struct DoubleKagome2 <: AbstractLattice
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

function DoubleKagome2(
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
    a2 = [-0.5a, 0.5√3a]

    # coordinates of each site in the unit cell
    r1 = [0.0, 0.0]
    r2 = 0.25 * a1
    r3 = 0.5 * [-a2[1], a2[2]]
    r4 = 0.5 * a1
    r5 = 0.75 * a1
    r6 = [2.5t, 0.5√3t]
    r = [r1, r2, r3, r4, r5, r6]

    return DoubleKagome2(n1, n2, t, a1, a2, r, PBC, antiPBC)
end

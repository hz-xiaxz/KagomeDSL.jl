using LinearAlgebra
abstract type AbstractLattice end

# note the license since I'm using code from BloqadeQMC

struct Kagome <: AbstractLattice
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
    distance_matrix::Array{Float64,2}
end

function create_distance_matrix(
    n1::Int,
    n2::Int,
    a1::Vector{Float64},
    a2::Vector{Float64},
    r::Array{Array{Float64,1},1},
    PBC::Tuple{Bool,Bool};
    trunc::Float64 = Inf,
)
    PBC1, PBC2 = PBC

    a2_length = sqrt(a2[1]^2 + a2[2]^2) # length of a2 vector
    θ = acos(a2[1] / a2_length) # a2 angle from horizontal

    N = n1 * n2 * length(r) # total number of sites
    dij = zeros(Float64, N, N) # distance matrix
    num_cells = n1 * n2 # number of repitions of the unit cell

    a2_length = sqrt(a2[1]^2 + a2[2]^2) # length of a2 vector
    θ = acos(a2[1] / a2_length) # a2 angle from horizontal

    # NOT ENDING ON num_cells-1 BECAUSE WE NEED INTER-CELL BONDS
    for i = 1:num_cells
        cellp1_x = rem(i, n1) > 0 ? rem(i, n1) - 1 : n1 - 1
        cellp1_y = rem(i, n1) > 0 ? div(i, n1) : div(i, n1) - 1

        cell1_x = cellp1_x * a1[1] + cellp1_y * a2_length * cos(θ)
        cell1_y = cellp1_y * a2_length * sin(θ)
        cell1 = [cell1_x, cell1_y]

        # NOT STARTING FROM i+1 BECAUSE WE NEED INTER-CELL BONDS
        for j = i:num_cells
            cellp2_x = rem(j, n1) > 0 ? rem(j, n1) - 1 : n1 - 1
            cellp2_y = rem(j, n1) > 0 ? div(j, n1) : div(j, n1) - 1

            cell2_x = cellp2_x * a1[1] + cellp2_y * a2_length * cos(θ)
            cell2_y = cellp2_y * a2_length * sin(θ)
            cell2 = [cell2_x, cell2_y]

            for site_i = 1:length(r)
                site_num_i = site_i + length(r) * (i - 1)
                ri = r[site_i] + cell1

                for site_j = 1:length(r)
                    site_num_j = site_j + length(r) * (j - 1)

                    # ensure that there are no diagonal entries
                    if site_num_i == site_num_j
                        continue
                    end

                    rj = r[site_j] + cell2

                    Δ = ri - rj
                    Δx, Δy = Δ[1], Δ[2]
                    d_np = sqrt(Δx^2 + Δy^2)

                    # checks for PBCs
                    if PBC1 & PBC2
                        # periodic in both lattice directions
                        # essentially, find smallest distance

                        # distances for periodicity in a1
                        rj .+= a1 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_1 = sqrt(Δx^2 + Δy^2)

                        rj .-= 2 * a1 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_2 = sqrt(Δx^2 + Δy^2)

                        # undo change
                        rj .+= a1 * n1

                        # distances for periodicity in a2
                        rj .+= a2 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_3 = sqrt(Δx^2 + Δy^2)

                        rj .-= 2 * a2 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_4 = sqrt(Δx^2 + Δy^2)

                        # undo change
                        rj .+= a2 * n1

                        # take the minimum distance
                        d = min(d_np, d_1, d_2, d_3, d_4)

                    elseif PBC1 & !PBC2
                        # distances for periodicity in a1
                        rj .+= a1 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_1 = sqrt(Δx^2 + Δy^2)

                        rj .-= 2 * a1 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_2 = sqrt(Δx^2 + Δy^2)

                        # undo change
                        rj .+= a1 * n1

                        d = min(d_np, d_1, d_2)

                    elseif PBC2 & !PBC1
                        # distances for periodicity in a2
                        rj .+= a2 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_3 = sqrt(Δx^2 + Δy^2)

                        rj .-= 2 * a2 * n1
                        Δ = ri - rj
                        Δx, Δy = Δ[1], Δ[2]
                        d_4 = sqrt(Δx^2 + Δy^2)

                        # undo change
                        rj .+= a2 * n1

                        # take the minimum distance
                        d = min(d_np, d_3, d_4)

                    else
                        d = d_np
                    end

                    dij[site_num_i, site_num_j] = d <= trunc ? d : 0.0

                end
            end
        end
    end
    return dij
end

create_distance_matrix(
    n1::Int,
    n2::Int,
    a1::Vector{Float64},
    a2::Vector{Float64},
    r::Array{Array{Float64,1},1},
    PBC::Bool;
    trunc::Float64 = Inf,
) = create_distance_matrix(n1, n2, a1, a2, r, (PBC, PBC); trunc = trunc)

function Kagome(t::Float64, n1::Int, n2::Int, PBC::Tuple{Bool,Bool}; trunc::Float64 = Inf)
    a = 2 * t
    a1 = [a, 0.0]
    a2 = [a * 0.5, a * sqrt(3) * 0.5]

    # coordinates of each site in the unit cell
    r1 = [0.0, 0.0]
    r2 = 0.5 * a1
    r3 = 0.5 * a2
    r = [r1, r2, r3]

    distance_matrix = create_distance_matrix(n1, n2, a1, a2, r, PBC; trunc = trunc)

    return Kagome(n1, n2, t, a1, a2, r, PBC, distance_matrix)
end

Kagome(t::Float64, n1::Int, n2::Int, PBC::Bool; trunc::Float64 = Inf) =
    Kagome(t, n1, n2, (PBC, PBC); trunc = trunc)

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
    distance_matrix::Array{Float64,2}
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
    trunc::Float64 = Inf,
)
    @assert n1 % 2 == 0
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

    distance_matrix = create_distance_matrix(n1 ÷ 2, n2, a1, a2, r, PBC; trunc = trunc)

    return DoubleKagome(n1, n2, t, a1, a2, r, PBC, distance_matrix)
end

DoubleKagome(t::Float64, n1::Int, n2::Int, PBC::Bool; trunc::Float64 = Inf) =
    DoubleKagome(t, n1, n2, (PBC, PBC); trunc = trunc)

function nearestNeighbor(lat::AbstractLattice)
    # select only the upper triangular part of the distance matrix
    distance_matrix = triu(lat.distance_matrix)
    # extract all the indices of the minimum distance elements
    indices = findall(x -> x ≈ lat.t, distance_matrix)
    return indices
end


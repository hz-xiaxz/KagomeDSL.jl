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
    χ::Float64,
    link::Dict,
)
    bond_num = bondNum(PBCs1, PBCs2)
    if haskey(link, bond_num)
        # a bit hacky `=`
        Tunneling[s1, s2] = χ * link[bond_num]
    end

    bond_num2 = bondNum(PBCs2, PBCs1)
    if haskey(link, bond_num2)
        Tunneling[s2, s1] = χ * link[bond_num2]
    end
end

"""
    Hmat(lat::AbstractLattice, χ::Float64, N_up::Int, N_down::Int)
"""
function Hmat(lat::T, χ::Float64, N_up::Int, N_down::Int) where {T<:AbstractLattice}

    n1 = lat.n1
    n2 = lat.n2
    ns = n1 * n2 * 3
    N = N_up + N_down

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
    # TODO: add back link here
    link_inter = Dict(
        (3, 3n1 + 1) => -1,
        (3, 3n1 - 1) => -1,
        (6, 3n1 + 2) => -1,
        (6, 3n1 + 4) => 1,
        (5, 7) => 1,
        (1, -1) => 1,
        (1,3-3n1) => -1,
        (2, 6-3n1) => -1,
        (4,6-3n1) => 1,
        (5,9-3n1) => -1,

    )

    for bond in nearestNeighbor(lat)
        s1, s2 = Tuple(bond)
        # consider in-cell cases
        if (s1 - 1) ÷ 6 == (s2 - 1) ÷ 6
            setTunnel!(tunneling, s1, s2, s1, s2, χ, link_in)
        end
        # consider inter-cell cases
        # consider PBC

        # start from OBC

        if (s1 - 1) ÷ 6 != (s2 - 1) ÷ 6
            if lat.PBC == (false, false)
                setTunnel!(tunneling, s1, s2, s1, s2, χ, link_inter)

            elseif lat.PBC == (true, false)
                for PBCs1 in (s1 - 3n1, s1 + 3n1, s1)
                    for PBCs2 in (s2 - 3n1, s2 + 3n1, s2)
                        setTunnel!(tunneling, s1, s2, PBCs1, PBCs2, χ, link_inter)
                    end
                end

            elseif lat.PBC == (false, true)
                for PBCs1 in (s1 - 3ns, s1 + 3ns, s1)
                    for PBCs2 in (s2 - 3ns, s2 + 3ns, s2)
                        setTunnel!(tunneling, s1, s2, PBCs1, PBCs2, χ, link_inter)
                    end
                end

            elseif lat.PBC == (true, true)
                for PBCs1 in (s1 - 3n1, s1 + 3n1, s1 - 3ns, s1 + 3ns, s1)
                    for PBCs2 in (s1 - 3n1, s1 + 3n1, s2 - 3ns, s2 + 3ns, s2)
                        setTunnel!(tunneling, s1, s2, PBCs1, PBCs2, χ, link_inter)
                    end
                end
            end
        end
    end
    return tunneling
end


function get_sublattice_indices(lat::DoubleKagome)
    ns = lat.n1 * lat.n2 * 3
    num_unit_cells = ns ÷ 6

    sublattices = [Int[] for _ = 1:6]

    for i = 0:num_unit_cells-1
        for j = 1:6
            push!(sublattices[j], 6 * i + j)
        end
    end

    return sublattices
end

function get_reciprocal_vectors(lat::DoubleKagome)
    a1 = lat.a1
    a2 = lat.a2

    # Using the formula for 2D reciprocal lattice vectors
    # b1 = 2π * [a2_y, -a2_x] / (a1_x*a2_y - a1_y*a2_x)
    # b2 = 2π * [-a1_y, a1_x] / (a1_x*a2_y - a1_y*a2_x)

    det = a1[1] * a2[2] - a1[2] * a2[1]
    b1 = 2π / det * [a2[2], -a2[1]]
    b2 = 2π / det * [-a1[2], a1[1]]

    return b1, b2
end

function get_K_points(lat::DoubleKagome)
    b1, b2 = get_reciprocal_vectors(lat)
    ideal_K1 = (2 / 3) * b1 + (1 / 3) * b2
    ideal_K2 = (1 / 3) * b1 + (2 / 3) * b2

    n1_eff = lat.n1 ÷ 2
    n2 = lat.n2

    offset1 = lat.antiPBC[1] ? 0.5 : 0.0
    offset2 = lat.antiPBC[2] ? 0.5 : 0.0

    min_dist1 = Inf
    min_dist2 = Inf
    closest_K1 = ideal_K1
    closest_K2 = ideal_K2

    for m = -n1_eff:n1_eff
        for k = -n2:n2
            q = ((m + offset1) / n1_eff) * b1 + ((k + offset2) / n2) * b2

            dist1 = norm(q - ideal_K1)
            if dist1 < min_dist1
                min_dist1 = dist1
                closest_K1 = q
            end

            dist2 = norm(q - ideal_K2)
            if dist2 < min_dist2
                min_dist2 = dist2
                closest_K2 = q
            end
        end
    end

    return closest_K1, closest_K2
end

function spin_structure_factor(
    sz_values::Vector{Float64},
    q::Vector{Float64},
    sites::Vector{Int},
    lat::DoubleKagome,
)
    factor = 0.0 + 0.0im
    for i in sites
        r = get_site_coord(lat, i)
        factor += sz_values[i] * exp(im * dot(q, r))
    end
    return abs2(factor) / length(sites)
end

function calculate_S_xy(mc::MC)
    ns = length(mc.kappa_up)
    S_xy = zeros(ComplexF64, ns, ns)

    for i = 1:ns
        for j = 1:ns
            if i == j
                # For a single spin-1/2 site, <S_i^x S_i^x + S_i^y S_i^y> = <S_i^2 - (S_i^z)^2>
                # = 3/4 - 1/4 = 1/2.
                # This assumes single occupancy, which is the case in a Mott insulator.
                S_xy[i, i] = 0.5
                continue
            end
            # Splus Sminus case
            if mc.kappa_down[i] != 0 && mc.kappa_up[j] != 0
                l_up = mc.kappa_up[j]
                K_up = i
                l_down = mc.kappa_down[i]
                K_down = j
                # think there is a minus sign
                S_xy[i, j] = -0.5 * mc.W_up[K_up, l_up] * mc.W_down[K_down, l_down]
                # Sminus Splus case
            elseif mc.kappa_up[i] != 0 && mc.kappa_down[j] != 0
                l_up = mc.kappa_up[i]
                K_up = j
                l_down = mc.kappa_down[j]
                K_down = i
                # think there is a minus sign
                S_xy[i, j] = -0.5 * mc.W_up[K_up, l_up] * mc.W_down[K_down, l_down]
            else
                S_xy[i, j] = 0.0
            end

            S_xy[i, j] = C_xy_ij
        end
    end
    return S_xy
end

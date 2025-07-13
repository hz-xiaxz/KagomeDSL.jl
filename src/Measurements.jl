
function get_sublattice_indices(lat::DoubleKagome)
    ns = lat.n1 * lat.n2 * 3
    num_unit_cells = ns ÷ 6
    
    sublattices = [Int[] for _ in 1:6]
    
    for i in 0:num_unit_cells-1
        for j in 1:6
            push!(sublattices[j], 6*i + j)
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
    
    det = a1[1]*a2[2] - a1[2]*a2[1]
    b1 = 2π/det * [a2[2], -a2[1]]
    b2 = 2π/det * [-a1[2], a1[1]]
    
    return b1, b2
end

function get_K_points(lat::DoubleKagome)
    b1, b2 = get_reciprocal_vectors(lat)
    ideal_K1 = (2/3)*b1 + (1/3)*b2
    ideal_K2 = (1/3)*b1 + (2/3)*b2

    n1_eff = lat.n1 ÷ 2
    n2 = lat.n2

    offset1 = lat.antiPBC[1] ? 0.5 : 0.0
    offset2 = lat.antiPBC[2] ? 0.5 : 0.0

    min_dist1 = Inf
    min_dist2 = Inf
    closest_K1 = ideal_K1
    closest_K2 = ideal_K2

    for m in -n1_eff:n1_eff
        for k in -n2:n2
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

function spin_structure_factor(sz_values::Vector{Float64}, q::Vector{Float64}, sites::Vector{Int}, lat::DoubleKagome)
    factor = 0.0 + 0.0im
    for i in sites
        r = get_site_coord(lat, i)
        factor += sz_values[i] * exp(im * dot(q, r))
    end
    return abs2(factor) / length(sites)
end

function calculate_S_xy(mc::MC, q::Vector{Float64})
    S_xy = 0.0 + 0.0im
    ns = length(mc.kappa_up)
    lat = mc.Ham.lat

    for i in 1:ns
        for j in 1:ns
            xprime = Dict{ConfigKey,Float64}()
            spinInteraction!(xprime, mc.kappa_up, mc.kappa_down, i, j)
            
            C_xy_ij = 0.0
            for (conf, coeff) in pairs(xprime)
                # The coefficient from spinInteraction! is -1/2.
                # Sx_i*Sx_j + Sy_i*Sy_j = (1/2) * (S+_i*S-_j + S-_i*S+_j)
                # So the value is -1 * coeff * ratio.
                update_up = mc.W_up[conf[1], conf[2]]
                update_down = mc.W_down[conf[3], conf[4]]
                ratio = update_up * update_down
                C_xy_ij += -coeff * ratio
            end

            if C_xy_ij != 0.0
                r_i = get_site_coord(lat, i)
                r_j = get_site_coord(lat, j)
                S_xy += C_xy_ij * exp(-im * dot(q, r_i - r_j))
            end
        end
    end
    return real(S_xy) / ns
end

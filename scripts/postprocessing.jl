using CairoMakie
using Carlo.ResultTools
using DataFrames
using Measurements
using KagomeDSL
using LinearAlgebra # For norm()

# --- Custom Kagome Lattice Definition (3-site unit cell) ---
abstract type AbstractLattice end

struct Kagome <: AbstractLattice
    n1::Int # number of repetitions in a1 direction
    n2::Int # number of repetitions in a2 direction
    t::Float64 # parameter that defines the equilateral triangle side length
    a1::Array{Float64,1}
    a2::Array{Float64,1}
    r::Array{Array{Float64,1},1} # coordinates of sites inside a unit cell
    PBC::Tuple{Bool,Bool}
    antiPBC::Tuple{Bool,Bool}
end

function Kagome(
    t::Float64,
    n1::Int,
    n2::Int,
    PBC::Tuple{Bool,Bool};
    antiPBC::Tuple{Bool,Bool} = (false, false),
)
    a1 = [2.0 * t, 0.0]
    a2 = [1.0 * t, sqrt(3.0) * t]

    r1 = [0.0, 0.0]
    r2 = [1.0 * t, 0.0]
    r3 = [0.5 * t, sqrt(3.0) / 2.0 * t]
    r = [r1, r2, r3]

    return Kagome(n1, n2, t, a1, a2, r, PBC, antiPBC)
end

ns(lat::AbstractLattice) = lat.n1 * lat.n2 * length(lat.r)
# --- End Custom Kagome Lattice Definition ---

# --- Your Existing Code (Slightly tidied) ---
# Load data
path = joinpath(@__DIR__, "../data/LL-4x4.results.json")
df = DataFrame(ResultTools.dataframe(path))
Sz = df[1, :Sz]
S_xy_values = df[1, :S_xy]

t = 1.0
n1 = 4
n2 = 4
# Note: KagomeDSL seems to create a 2*n1 x n2 lattice.
# For a 4x4 unit cell Kagome, n1=4, n2=4 gives an 8x4 system.
# Let's assume the lattice generation is correct for your system.
lattice = Kagome(t, n1, n2, (true, true))
sites = [
    (i - 1) * lattice.a1 + (j - 1) * lattice.a2 + r_ for j = 1:lattice.n2 for
    i = 1:lattice.n1 for r_ in lattice.r
]
x_coords = [s[1] for s in sites]
y_coords = [s[2] for s in sites]
num_sites = length(sites)

# Extract numerical values
sz_values = Measurements.value.(Sz)
s_xy_values = Measurements.value.(S_xy)

# --- Create the Figure with 2 Panels ---
fig = Figure()



# --- NEW: Panel 2: Inferred Spin Directions (Vector Plot) ---
ax = Axis(
    fig[1, 1],
    aspect = DataAspect(),
    title = "Inferred Spin Directions (Pinned at Red Site)",
)

annotation!(
    ax,
    x_coords,
    y_coords,
    text = string.(1:length(sites)),
    color = :black,
    fontsize = 8,
)

# CHOOSE YOUR REFERENCE SITE HERE! A site near the center is usually best.
# ref_site_idx = 1 # Let's pick site 25 as an example

# --- MODIFICATION: Account for non-zero Sz in angle calculation ---
# We want to calculate cos(Δθ) = <S_ref_xy . S_j_xy> / (|S_ref_xy| |S_j_xy|)
# We approximate |S_i_xy| using <S_i_z>.
# For a spin S, |S_i_xy|^2 = S(S+1) - S_iz^2.
# We approximate <|S_i_xy|^2> ≈ S(S+1) - <S_iz^2>
# And we approximate <S_iz^2> ≈ <S_iz>^2 + Var(S_iz)
S = 0.5 # Assume S=1/2

# CHOOSE YOUR REFERENCE SITE HERE! A site near the center is usually best.
ref_site_idx = 1# Let's pick site 25 as an example

# The auto-correlation of the reference site is used for scaling arrow lengths
C_ref_auto = real(s_xy_values[ref_site_idx, ref_site_idx])
# --- END MODIFICATION ---

# Calculate the vector components (u, v) for each site
u = zeros(num_sites)
v = zeros(num_sites)

# Function to get Kagome neighbors based on lattice structure
function calculate_intra_unitcell_angles(
    current_site_idx::Int,
    known_angle::Float64,
    s_xy_values::Matrix{Float64},
    lattice::Kagome,
    num_sites::Int,
)
    intra_unitcell_angles = Dict{Int,Float64}()

    # Convert linear site_idx to (i, j, r_idx) (0-based for calculations)
    r_idx_0based = (current_site_idx - 1) % length(lattice.r)
    temp = (current_site_idx - 1) ÷ length(lattice.r)
    i_0based = temp % lattice.n1
    j_0based = temp ÷ lattice.n1

    # Identify other r_idx within the same unit cell
    other_r_idx_0based = filter(x -> x != r_idx_0based, 0:length(lattice.r)-1)

    for target_r_idx_0based in other_r_idx_0based
        # Calculate linear index of the intra-unit-cell neighbor
        neighbor_linear_idx =
            j_0based * lattice.n1 * length(lattice.r) +
            i_0based * length(lattice.r) +
            target_r_idx_0based +
            1

        if 1 <= neighbor_linear_idx <= num_sites
            # Get correlation with the current site
            corr_val = real(s_xy_values[current_site_idx, neighbor_linear_idx])

            # Calculate delta_theta
            delta_theta = acos(corr_val / 0.25)

            # Apply sign flip if it's the "third" site in its unit cell (r3)
            if target_r_idx_0based + 1 == 3 # Convert back to 1-based for comparison with r[3]
                delta_theta = -delta_theta
            end

            # Calculate the absolute angle of the neighbor
            neighbor_angle = known_angle + delta_theta

            intra_unitcell_angles[neighbor_linear_idx] = neighbor_angle
        end
    end

    return intra_unitcell_angles
end

# Function to generate all nearest-neighbor bonds for the entire lattice
function generate_all_nn_bonds(lattice::Kagome, num_sites::Int)
    all_nn_bonds = Set{Tuple{Int,Int}}()

    # Helper to convert (i, j, r_idx) to linear index, handling OBC
    function to_linear_idx(i, j, r_idx)
        if !(0 <= i < lattice.n1) || !(0 <= j < lattice.n2)
            return -1 # Out of bounds
        end
        return j * lattice.n1 * length(lattice.r) + i * length(lattice.r) + r_idx + 1
    end

    for site_idx = 1:num_sites
        # Convert linear site_idx to (i, j, r_idx) (0-based for calculations)
        r_idx_0based = (site_idx - 1) % length(lattice.r)
        i_0based = ((site_idx - 1) ÷ length(lattice.r)) % lattice.n1
        j_0based = ((site_idx - 1) ÷ length(lattice.r)) ÷ lattice.n1

        # Define neighbor connections based on r_idx_0based
        if r_idx_0based == 0 # r[1] site
            # Within same unit cell
            potential_neighbors = [
                to_linear_idx(i_0based, j_0based, 1), # r[2]
                to_linear_idx(i_0based, j_0based, 2), # r[3]
                # Inter-unit cell
                to_linear_idx(i_0based - 1, j_0based, 1), # r[2] of (i-1, j) cell
                to_linear_idx(i_0based, j_0based - 1, 2), # r[3] of (i, j-1) cell
            ]
        elseif r_idx_0based == 1 # r[2] site
            # Within same unit cell
            potential_neighbors = [
                to_linear_idx(i_0based, j_0based, 0), # r[1]
                to_linear_idx(i_0based, j_0based, 2), # r[3]
                # Inter-unit cell
                to_linear_idx(i_0based + 1, j_0based, 0), # r[1] of (i+1, j) cell
                to_linear_idx(i_0based + 1, j_0based - 1, 2), # r[3] of (i+1, j-1) cell
            ]
        elseif r_idx_0based == 2 # r[3] site
            # Within same unit cell
            potential_neighbors = [
                to_linear_idx(i_0based, j_0based, 0), # r[1]
                to_linear_idx(i_0based, j_0based, 1), # r[2]
                # Inter-unit cell
                to_linear_idx(i_0based, j_0based + 1, 0), # r[1] of (i, j+1) cell
                to_linear_idx(i_0based - 1, j_0based + 1, 1), # r[2] of (i-1, j+1) cell
            ]
        else
            potential_neighbors = Int[] # Should not happen
        end

        for neighbor_linear_idx in potential_neighbors
            if 1 <= neighbor_linear_idx <= num_sites
                # Ensure bond is added only once (smaller index first)
                if site_idx < neighbor_linear_idx
                    push!(all_nn_bonds, (site_idx, neighbor_linear_idx))
                else
                    push!(all_nn_bonds, (neighbor_linear_idx, site_idx))
                end
            end
        end
    end
    return sort(collect(all_nn_bonds))
end

# Function to calculate angles of other sites within the same unit cell
function calculate_intra_unitcell_angles(
    current_site_idx::Int,
    known_angle::Float64,
    s_xy_values::AbstractMatrix,
    lattice::Kagome,
    num_sites::Int,
)
    intra_unitcell_angles = Dict{Int,Float64}()

    # Convert linear site_idx to (i, j, r_idx) (0-based for calculations)
    r_idx_0based = (current_site_idx - 1) % length(lattice.r)
    temp = (current_site_idx - 1) ÷ length(lattice.r)
    i_0based = temp % lattice.n1
    j_0based = temp ÷ lattice.n1

    # Identify other r_idx within the same unit cell
    other_r_idx_0based = filter(x -> x != r_idx_0based, 0:length(lattice.r)-1)

    for target_r_idx_0based in other_r_idx_0based
        # Calculate linear index of the intra-unit-cell neighbor
        neighbor_linear_idx =
            j_0based * lattice.n1 * length(lattice.r) +
            i_0based * length(lattice.r) +
            target_r_idx_0based +
            1

        if 1 <= neighbor_linear_idx <= num_sites
            # Get correlation with the current site
            corr_val = real(s_xy_values[current_site_idx, neighbor_linear_idx])
            # Calculate delta_theta
            delta_theta = acos(corr_val / 0.25)

            # Apply sign flip if it's the "third" site in its unit cell (r3)
            if target_r_idx_0based + 1 == 3 # Convert back to 1-based for comparison with r[3]
                delta_theta = -delta_theta
            end

            # Calculate the absolute angle of the neighbor
            neighbor_angle = known_angle + delta_theta

            intra_unitcell_angles[neighbor_linear_idx] = neighbor_angle
        end
    end

    return intra_unitcell_angles
end

# Initialize for iterative angle propagation
known_angles = Dict{Int,Float64}()

# Pin the initial reference site (site 1), which is an r[1] site
ref_site_idx = 1
theta_ref_initial = atan(-1 / 2, -sqrt(3) / 2) # Initial angle for site 1
known_angles[ref_site_idx] = theta_ref_initial

# Generate all nearest-neighbor bonds once
all_nn_bonds = generate_all_nn_bonds(lattice, num_sites)

# Iteratively propagate angles
let changed = true
    while changed
        changed = false
        for (s1, s2) in all_nn_bonds
            if s1 > s2
                continue # Ensure we process each bond only once
            end
            # Case 1: s1's angle is known, s2's is not
            if haskey(known_angles, s1) && !haskey(known_angles, s2)
                current_site_idx = s1
                neighbor_site_idx = s2
                current_site_angle = known_angles[current_site_idx]

                corr_val = real(s_xy_values[current_site_idx, neighbor_site_idx])
                delta_theta = acos(corr_val / 0.25)

                r_idx_1based_neighbor = ((neighbor_site_idx - 1) % 3) + 1
                if r_idx_1based_neighbor == 3
                    delta_theta = -delta_theta
                end

                known_angles[neighbor_site_idx] = current_site_angle + delta_theta
                changed = true

                # If the newly determined site is an r[1] site, calculate its intra-unit-cell partners
                if ((neighbor_site_idx - 1) % 3) + 1 == 1
                    intra_angles = calculate_intra_unitcell_angles(
                        neighbor_site_idx,
                        known_angles[neighbor_site_idx],
                        s_xy_values,
                        lattice,
                        num_sites,
                    )
                    for (idx, angle) in intra_angles
                        if !haskey(known_angles, idx)
                            known_angles[idx] = angle
                            changed = true
                        end
                    end
                end
            end
        end
    end
end

# Populate u and v arrays from known_angles
for j = 1:num_sites
    if haskey(known_angles, j)
        unitcell = (j - 1) ÷ 3 + 1
        if unitcell % n1 == 2
            # Calculate the unit vector components based on the known angle
            known_angles[j] -= 4π / 3
        elseif unitcell % n1 == 3
            # Calculate the unit vector components based on the known angle
            known_angles[j] -= 2π / 3
        end
        u[j] = cos(known_angles[j])
        v[j] = sin(known_angles[j])
    else
        # This case should ideally not be reached if the lattice is connected
        # and all sites are reachable from site 1.
        # If it is reached, it means some sites were not processed.
        # For now, we set them to zero and issue a warning.
        u[j] = 0.0
        v[j] = 0.0
        @warn "Site $j was not processed. Check lattice connectivity or BFS logic."
    end
end


# Plot the bonds
for (s1, s2) in all_nn_bonds
    lines!(
        ax,
        [x_coords[s1], x_coords[s2]],
        [y_coords[s1], y_coords[s2]],
        color = :gray,
        linewidth = 0.5,
    )
end

# Plot the vector field
arrows2d!(
    ax,
    x_coords,
    y_coords,
    u,
    v,
    tiplength = 10,
    normalize = true,
    lengthscale = 0.5, # Adjust this to make arrows larger/smaller
    tipcolor = :blue,
    shaftcolor = :blue,
)
# Highlight the reference spin in red
arrows2d!(
    ax,
    [x_coords[ref_site_idx]],
    [y_coords[ref_site_idx]],
    [u[ref_site_idx]],
    [v[ref_site_idx]],
    tiplength = 15,
    lengthscale = 0.5,
    tipcolor = :red,
    shaftcolor = :red,
)

scatter!(ax, x_coords, y_coords, color = :lightgray, markersize = 5) # Show site positions faintly

hidedecorations!(ax)
hidespines!(ax)


# Display the final figure
save("./correlation.png", fig)
fig

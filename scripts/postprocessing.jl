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
Sz = df[5, :Sz]
S_xy_values = df[5, :S_xy]

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
num_sites = 3 * n1 * n2

# Initialize for iterative angle propagation
known_angles = zeros(Float64, num_sites)

function cal_unitcell(
    first_site_angle::Float64,
    first_site_idx::Int,
    s_xy_values::AbstractMatrix,
)
    corr12 = real(s_xy_values[first_site_idx, first_site_idx+1])
    corr13 = real(s_xy_values[first_site_idx, first_site_idx+2])
    angle12 = acos(corr12 / 0.25) # ∈ [0, π]
    angle13 = acos(corr13 / 0.25) # ∈ [0, π]
    angle2 = first_site_angle + angle12
    angle3 = first_site_angle - angle13
    return angle2, angle3
end

first_site_angles = zeros(Float64, n1 * n2)

for i = 1:(n1*n2)
    first_site = (i - 1) * 3 + 1
    if i == 1
        first_site_angles[i] = atan(-1 / 2, -sqrt(3) / 2) # Initial angle for site 1
        angle2, angle3 = cal_unitcell(first_site_angles[i], first_site, s_xy_values)
        known_angles[first_site] = first_site_angles[i]
        known_angles[first_site+1] = angle2 
        known_angles[first_site+2] = angle3 
        continue
    end
    # deal with first site
    # deal with every beginning of the row
    if (i - 1) % n1 == 0
        linked_site = first_site - 3 * n1 + 2
        corr = real(s_xy_values[linked_site, first_site])
        anglediff = acos(corr / 0.25)
        first_site_angles[i] = known_angles[linked_site] + anglediff
    else
        linked_site = first_site - 2
        corr = real(s_xy_values[linked_site, first_site])
        anglediff = acos(corr / 0.25)
        first_site_angles[i] = known_angles[linked_site] - anglediff
    end

    angle2, angle3 = cal_unitcell(first_site_angles[i], first_site, s_xy_values)
    known_angles[first_site] = first_site_angles[i]
    known_angles[first_site+1] = angle2 
    known_angles[first_site+2] = angle3 
end
u = zeros(Float64, num_sites)
v = zeros(Float64, num_sites)
for i = 1:num_sites
    u[i] = cos(known_angles[i])
    v[i] = sin(known_angles[i])
end
@show known_angles./π

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

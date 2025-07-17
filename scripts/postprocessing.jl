using CairoMakie
using Carlo.ResultTools
using DataFrames
using Measurements
using KagomeDSL
using LinearAlgebra # For norm()

# Load data
path = joinpath(@__DIR__, "../data/LL-4x4.results.json")
df = DataFrame(ResultTools.dataframe(path))
Sz = df[1, :Sz]
S_xy_values = df[1, :S_xy]

t = 1.0
n1 = 4
n2 = 4
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
s_xy_values = Measurements.value.(S_xy_values)

# --- Create the Figure with 2 Panels ---
fig = Figure()

ax = Axis(
    fig[1, 1],
    aspect = DataAspect(),
    title = "Inferred Spin Directions (Pinned at Red Site)",
)

annotation!(
    ax,
    x_coords,
    y_coords,
    text = string.(1:num_sites),
    color = :black,
    fontsize = 8,
)

# Calculate the vector components (u, v) for each site

u = zeros(num_sites)
v = zeros(num_sites)


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

for i = 1:num_sites
    u[i] = cos(known_angles[i])
    v[i] = sin(known_angles[i])
end
@show known_angles ./ π

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

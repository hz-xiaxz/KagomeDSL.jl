using CairoMakie
using Carlo.ResultTools
using DataFrames
using Measurements
using KagomeDSL
using LinearAlgebra # For norm()

# --- Your Existing Code (Slightly tidied) ---
# Load data
path = joinpath(@__DIR__, "../data/LL-4x4.results.json")
df = DataFrame(ResultTools.dataframe(path))
Sz = df[1, :Sz]
S_xy = df[1, :S_xy]

t = 1.0
n1 = 4
n2 = 4
# Note: KagomeDSL seems to create a 2*n1 x n2 lattice.
# For a 4x4 unit cell Kagome, n1=4, n2=4 gives an 8x4 system.
# Let's assume the lattice generation is correct for your system.
lattice = KagomeDSL.DoubleKagome(t, n1, n2, (true, true))
sites = [
    (i - 1) * lattice.a1 + (j - 1) * lattice.a2 + r_ for j = 1:lattice.n2 for
    i = 1:lattice.n1÷2 for r_ in lattice.r
]
x_coords = [s[1] for s in sites]
y_coords = [s[2] for s in sites]
num_sites = length(sites)

# Extract numerical values
sz_values = Measurements.value.(Sz)
s_xy_values = Measurements.value.(S_xy)

# --- Create the Figure with 3 Panels ---
fig = Figure()

# Panel 1: Sz expectation values (your original plot)
ax1 = Axis(fig[1, 1], aspect = DataAspect(), title = "⟨Sz⟩ Expectation Values")
scatter!(ax1, x_coords, y_coords, color = sz_values, colormap = :viridis)
annotation!(
    ax1,
    x_coords,
    y_coords,
    text = string.(1:length(sites)),
    color = :black,
    fontsize = 8,
)
Colorbar(fig[1, 1][1, 2], limits = extrema(sz_values), colormap = :viridis)
# Hiding decorations for a cleaner look
hidedecorations!(ax1)
hidespines!(ax1)


# Panel 2: S_xy correlations as bonds (your original plot)
ax2 = Axis(fig[1, 2], aspect = DataAspect(), title = "S_xy Correlations (Bond Plot)")
scatter!(ax2, x_coords, y_coords, color = :gray)
annotation!(
    ax2,
    x_coords,
    y_coords,
    text = string.(1:length(sites)),
    color = :black,
    fontsize = 8,
)
s_xy_max_abs = maximum(abs.(s_xy_values))
for i = 1:num_sites, j = (i+1):num_sites
    corr_val = s_xy_values[i, j]
    if abs(corr_val) > 0.01 # Threshold to avoid clutter
        lines!(
            ax2,
            [x_coords[i], x_coords[j]],
            [y_coords[i], y_coords[j]],
            color = corr_val,
            colormap = :viridis,
            colorrange = (-s_xy_max_abs, s_xy_max_abs),
            linewidth = 8 * abs(corr_val) / s_xy_max_abs,
        )
    end
end
Colorbar(
    fig[1, 2][1, 2],
    limits = (-s_xy_max_abs, s_xy_max_abs),
    colormap = :viridis,
    label = "⟨SᵢˣSⱼˣ+SᵢʸSⱼʸ⟩",
)
hidedecorations!(ax2)
hidespines!(ax2)


# --- NEW: Panel 3: Inferred Spin Directions (Vector Plot) ---
ax3 = Axis(
    fig[1, 3],
    aspect = DataAspect(),
    title = "Inferred Spin Directions (Pinned at Red Site)",
)

annotation!(
    ax3,
    x_coords,
    y_coords,
    text = string.(1:length(sites)),
    color = :black,
    fontsize = 8,
)

# CHOOSE YOUR REFERENCE SITE HERE! A site near the center is usually best.
ref_site_idx = 1 # Let's pick site 25 as an example

# The auto-correlation C(i,i) gives the theoretical maximum value, S(S+1)/2 or S^2
# Using it makes the normalization robust. For S=1/2, this is ~0.25
C_max = s_xy_values[ref_site_idx, ref_site_idx]
if C_max < 1e-6
    @warn "Auto-correlation at reference site is near zero. Check your data or ref_site_idx."
    C_max = 0.25# Fallback for S=1/2
end

# Calculate the vector components (u, v) for each site
u = zeros(num_sites)
v = zeros(num_sites)

for j = 1:num_sites
    if j == ref_site_idx
        # Pin the reference spin to point right
        u[j] = -1 / 2
        v[j] = -√3 / 2
    else
        # Get correlation with the reference site
        corr_val = s_xy_values[ref_site_idx, j]

        # Normalize correlation to get cos(Δθ)
        # clamp ensures the value is in [-1, 1] to avoid acos domain errors
        c_norm = clamp(corr_val / C_max, -1.0, 1.0)

        # The angle of spin j relative to the reference spin (which is at 0 rad)
        delta_theta = acos(c_norm)

        # The length of the arrow is proportional to the correlation strength
        arrow_length = abs(corr_val / C_max)

        # Convert polar (length, angle) to Cartesian (u, v) components
        u[j] = arrow_length * cos(delta_theta)
        v[j] = arrow_length * sin(delta_theta)
    end
end

# Plot the vector field
arrows2d!(
    ax3,
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
    ax3,
    [x_coords[ref_site_idx]],
    [y_coords[ref_site_idx]],
    [u[ref_site_idx]],
    [v[ref_site_idx]],
    tiplength = 15,
    lengthscale = 0.5,
    tipcolor = :red,
    shaftcolor = :red,
)

scatter!(ax3, x_coords, y_coords, color = :lightgray, markersize = 5) # Show site positions faintly

hidedecorations!(ax3)
hidespines!(ax3)


# Display the final figure
save("./correlation.png", fig)
fig

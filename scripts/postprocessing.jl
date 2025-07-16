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
Sz = df[6, :Sz]
S_xy = df[6, :S_xy]

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

# --- Create the Figure with 2 Panels ---
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


# --- NEW: Panel 2: Inferred Spin Directions (Vector Plot) ---
ax2 = Axis(
    fig[1, 2],
    aspect = DataAspect(),
    title = "Inferred Spin Directions (Pinned at Red Site)",
)

annotation!(
    ax2,
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
S_sq_plus_1 = S * (S + 1)

sz_exp_vals = Measurements.value.(Sz)
sz_uncertainties = Measurements.uncertainty.(Sz)
sz_sq_exp = sz_exp_vals .^ 2 + sz_uncertainties .^ 2

# Estimate of <|S_i_xy|^2> for each site
s_perp_sq = max.(0, S_sq_plus_1 .- sz_sq_exp)
# Estimate of sqrt(<|S_i_xy|^2>) for each site
s_perp = sqrt.(s_perp_sq)

# CHOOSE YOUR REFERENCE SITE HERE! A site near the center is usually best.
ref_site_idx = 1 # Let's pick site 25 as an example

# The auto-correlation of the reference site is used for scaling arrow lengths
C_ref_auto = s_xy_values[ref_site_idx, ref_site_idx]
if C_ref_auto < 1e-6
    @warn "Auto-correlation at reference site is near zero. Check your data or ref_site_idx."
    C_ref_auto = S_sq_plus_1 # Fallback for S=1/2, S(S+1)
end
# --- END MODIFICATION ---

# Calculate the vector components (u, v) for each site
u = zeros(num_sites)
v = zeros(num_sites)

# The reference spin is pinned at angle theta_ref. We add the relative
# angle delta_theta to get the absolute angle of the current spin.
# Note: There is an inherent ambiguity in the sign of delta_theta.
# We consistently choose the positive root from acos.
theta_ref = atan(-√3 / 2, -1 / 2)

for j = 1:num_sites
    if j == ref_site_idx
        # Pin the reference spin
        u[j] = -1 / 2
        v[j] = -√3 / 2
    else
        # Get correlation with the reference site
        corr_val = s_xy_values[ref_site_idx, j]

        # --- MODIFICATION: Use new normalization for cos(Δθ) ---
        # Normalize correlation to get cos(Δθ)
        # The normalization factor is |S_ref_xy| * |S_j_xy|
        norm_factor = s_perp[ref_site_idx] * s_perp[j]

        c_norm = if norm_factor < 1e-6
            0.0
        else
            # clamp ensures the value is in [-1, 1] to avoid acos domain errors
            clamp(corr_val / norm_factor, -1.0, 1.0)
        end
        # --- END MODIFICATION ---

        # The angle of spin j relative to the reference spin
        delta_theta = acos(c_norm)

        # The length of the arrow is proportional to the correlation strength
        # Let's normalize by the auto-correlation of the reference site
        arrow_length = abs(corr_val / C_ref_auto)

        # Add the relative angle to the reference angle to get the absolute angle
        theta_j = theta_ref + delta_theta

        # Convert polar (length, angle) to Cartesian (u, v) components
        u[j] = arrow_length * cos(theta_j)
        v[j] = arrow_length * sin(theta_j)
    end
end

# Plot the vector field
arrows2d!(
    ax2,
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
    ax2,
    [x_coords[ref_site_idx]],
    [y_coords[ref_site_idx]],
    [u[ref_site_idx]],
    [v[ref_site_idx]],
    tiplength = 15,
    lengthscale = 0.5,
    tipcolor = :red,
    shaftcolor = :red,
)

scatter!(ax2, x_coords, y_coords, color = :lightgray, markersize = 5) # Show site positions faintly

hidedecorations!(ax2)
hidespines!(ax2)


# Display the final figure
save("./correlation.png", fig)
fig

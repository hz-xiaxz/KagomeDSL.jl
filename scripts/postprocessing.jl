using CairoMakie
using Carlo.ResultTools
using DataFrames
using Measurements
using KagomeDSL
using LinearAlgebra # For norm()
using MakiePublication
set_theme!(theme_web())

# Load data
path = joinpath(@__DIR__, "../data/LL-4x4.results.json")
df = DataFrame(ResultTools.dataframe(path))
Sz = Observable(df[1, :Sz])
S_xy_values = Observable(df[1, :S_xy])
imbalance = Observable(df[1, :imbalance])

t = 1.0
n1 = 4
n2 = 4
# For a 4x4 unit cell Kagome, n1=4, n2=4 gives an 8x4 system.
# Let's assume the lattice generation is correct for your system.
lattice = DoubleKagome(t, n1, n2, (true, true); antiPBC = (false, false))
sites = [
    (i - 1) * lattice.a1 + (j - 1) * lattice.a2 + r_ for j = 1:lattice.n2 for
    i = 1:(lattice.n1÷2) for r_ in lattice.r
]
x_coords = [s[1] for s in sites]
y_coords = [s[2] for s in sites]
num_sites = length(sites)


# Extract numerical values
sz_values = @lift(Measurements.value.($(Sz)))
s_xy_values = @lift(Measurements.value.($(S_xy_values)))

# --- Create the Figure with 2 Panels ---
fig = annotation(
    x_coords,
    y_coords,
    text = string.(1:num_sites),
    color = :black,
    fontsize = 8,
    axis = (
        title = @lift(
            "Inferred Spin Directions (Pinned at Red Site) imbalance:$($imbalance)"
        ),
    ),
)

# Calculate the vector components (u, v) for each site

u = Observable(zeros(num_sites))
v = Observable(zeros(num_sites))

function update_uv(num_sites, s_xy_values)
    u = zeros(Float64, num_sites)
    v = zeros(Float64, num_sites)

    # Initialize for iterative angle propagation
    known_angles = zeros(Float64, num_sites)
    s_xy_values = Measurements.value.(real.(s_xy_values))

    function cal_unitcell(first_site_angle::Float64, first_site_idx::Int, s_xy_values)
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
    return (u, v)
end

# Plot the bonds
for (s1, s2) in all_nn_bonds
    lines!(
        [x_coords[s1], x_coords[s2]],
        [y_coords[s1], y_coords[s2]],
        color = :gray,
        linewidth = 0.5,
    )
end

# Plot the vector field
arrows2d!(
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
    [x_coords[ref_site_idx]],
    [y_coords[ref_site_idx]],
    [-√3 / 2],
    [-1 / 2],
    tiplength = 15,
    lengthscale = 0.5,
    normalize = true,
    tipcolor = :red,
    shaftcolor = :red,
)

scatter!(x_coords, y_coords, color = :lightgray, markersize = 5) # Show site positions faintly

hidedecorations!()
hidespines!()

framerate = 1
timestamps = 1:10
record(fig, "XY.mp4", timestamps; framerate = framerate) do i
    # Update the observable values
    Sz[] = df[i, :Sz]
    S_xy_values[] = df[i, :S_xy]

    imbalance[] = df[i, :imbalance]
    u[], v[] = update_uv(num_sites, Measurements.values.(df[i, :S_xy]))
end

for i = 1:10
    Sz[] = df[i, :Sz]
    S_xy_values[] = df[i, :S_xy]
    imbalance[] = df[i, :imbalance]
    u[], v[] = update_uv(num_sites, Measurements.values.(df[i, :S_xy]))
    # Update the plot
    save(joinpath(@__DIR__, "../data/LL-4x4-imbalance-$i.png"),fig)
end
using Carlo.ResultTools
using DataFrames
using Measurements
using KagomeDSL
using LinearAlgebra
using Statistics
using CairoMakie
using Printf

# Load data
path = joinpath(@__DIR__, "../data/LL-4x4.results.json")
df = DataFrame(ResultTools.dataframe(path))

# Parameters
t = 1.0
n1 = 4
n2 = 4

# Generate lattice and site coordinates
lattice = DoubleKagome(t, n1, n2, (true, true); antiPBC = (false, false))
sites = [
    (i - 1) * lattice.a1 + (j - 1) * lattice.a2 + r_ for j = 1:lattice.n2 for
    i = 1:(lattice.n1÷2) for r_ in lattice.r
]
x_coords = [s[1] for s in sites]
y_coords = [s[2] for s in sites]
num_sites = length(sites)

println("Number of sites: $num_sites")
println("Lattice parameters: n1=$n1, n2=$n2")

# Function to calculate distance between two sites
function site_distance(i, j, x_coords, y_coords)
    return sqrt((x_coords[i] - x_coords[j])^2 + (y_coords[i] - y_coords[j])^2)
end

# Analyze S_xy vs distance for each measurement
function analyze_sxy_vs_distance(s_xy_matrix, x_coords, y_coords)
    distances = Float64[]
    correlations = Float64[]
    correlation_errors = Float64[]

    num_sites = length(x_coords)

    for i = 1:num_sites
        for j = (i+1):num_sites  # Only consider unique pairs
            dist = site_distance(i, j, x_coords, y_coords)
            corr_val = s_xy_matrix[i, j]

            push!(distances, dist)
            push!(correlations, real(Measurements.value(corr_val)))
            push!(correlation_errors, Measurements.uncertainty(corr_val))
        end
    end

    return distances, correlations, correlation_errors
end

# Function to bin data by distance
function bin_by_distance(distances, correlations, correlation_errors, num_bins = 20)
    min_dist = minimum(distances)
    max_dist = maximum(distances)
    bin_width = (max_dist - min_dist) / num_bins

    bin_centers = Float64[]
    bin_means = Float64[]
    bin_stds = Float64[]
    bin_counts = Int[]

    for i = 1:num_bins
        bin_start = min_dist + (i - 1) * bin_width
        bin_end = min_dist + i * bin_width
        bin_center = (bin_start + bin_end) / 2

        # Find points in this bin
        in_bin = (distances .>= bin_start) .& (distances .< bin_end)
        if i == num_bins  # Include the last point in the final bin
            in_bin = in_bin .| (distances .== max_dist)
        end

        if sum(in_bin) > 0
            bin_corrs = correlations[in_bin]
            push!(bin_centers, bin_center)
            push!(bin_means, mean(bin_corrs))
            push!(bin_stds, std(bin_corrs))
            push!(bin_counts, sum(in_bin))
        end
    end

    return bin_centers, bin_means, bin_stds, bin_counts
end

# Analyze first measurement
println("\nAnalyzing S_xy vs distance for first measurement...")
s_xy_matrix = df[1, :S_xy]
distances, correlations, correlation_errors =
    analyze_sxy_vs_distance(s_xy_matrix, x_coords, y_coords)

println("Total number of site pairs: $(length(distances))")
println("Distance range: $(minimum(distances)) to $(maximum(distances))")
println("Correlation range: $(minimum(correlations)) to $(maximum(correlations))")

# Bin the data
bin_centers, bin_means, bin_stds, bin_counts =
    bin_by_distance(distances, correlations, correlation_errors)

# Create plots
fig = Figure()

# First subplot: scatter plot
ax1 = Axis(
    fig[1, 1],
    xlabel = "Distance",
    ylabel = "S_xy correlation",
    title = "S_xy vs Distance (All Pairs)",
)
scatter!(ax1, distances, correlations, markersize = 4, alpha = 0.6)

# Second subplot: binned plot with error bars
ax2 = Axis(
    fig[2, 1],
    xlabel = "Distance",
    ylabel = "Average S_xy correlation",
    title = "Binned S_xy vs Distance",
)
errorbars!(ax2, bin_centers, bin_means, bin_stds, alpha = 0.7)
lines!(ax2, bin_centers, bin_means, linewidth = 2, alpha = 0.7)

# Save the plot
save(joinpath(@__DIR__, "../data/sxy_distance_analysis.png"), fig)
println("Plot saved to: ../data/sxy_distance_analysis.png")

# Print binned statistics
println("\nBinned analysis:")
println("Distance Range | Mean S_xy | Std S_xy | Count")
println("---------------|-----------|----------|------")
for i = 1:length(bin_centers)
    @printf(
        "%.3f - %.3f   | %8.5f | %8.5f | %5d\n",
        bin_centers[i] - (bin_centers[2] - bin_centers[1]) / 2,
        bin_centers[i] + (bin_centers[2] - bin_centers[1]) / 2,
        bin_means[i],
        bin_stds[i],
        bin_counts[i]
    )
end

# Analyze correlation decay
println("\nCorrelation decay analysis:")
near_neighbors = distances .< 2.0  # Adjust threshold as needed
far_neighbors = distances .> 4.0   # Adjust threshold as needed

if sum(near_neighbors) > 0 && sum(far_neighbors) > 0
    near_corr = mean(correlations[near_neighbors])
    far_corr = mean(correlations[far_neighbors])

    println("Near neighbors (< 2.0): Mean S_xy = $(near_corr)")
    println("Far neighbors (> 4.0): Mean S_xy = $(far_corr)")
    println("Decay ratio: $(far_corr/near_corr)")
end

# Save numerical results
results_data = DataFrame(
    distance = distances,
    sxy_correlation = correlations,
    sxy_error = correlation_errors,
)

using CSV
CSV.write(joinpath(@__DIR__, "../data/sxy_distance_data.csv"), results_data)
println("\nNumerical data saved to: ../data/sxy_distance_data.csv")

# Analyze time evolution if multiple measurements available
if nrow(df) > 1
    println("\nAnalyzing time evolution...")

    # Take a few representative measurements
    time_indices = min(10, nrow(df))
    time_fig = Figure()
    time_ax = Axis(
        time_fig[1, 1],
        title = "S_xy vs Distance Over Time",
        xlabel = "Distance",
        ylabel = "S_xy correlation",
    )

    for i = 1:time_indices
        s_xy_i = df[i, :S_xy]
        dist_i, corr_i, err_i = analyze_sxy_vs_distance(s_xy_i, x_coords, y_coords)
        bin_centers_i, bin_means_i, _, _ = bin_by_distance(dist_i, corr_i, err_i)

        lines!(
            time_ax,
            bin_centers_i,
            bin_means_i,
            label = "t=$i",
            alpha = 0.7,
            linewidth = 2,
        )
    end

    axislegend(time_ax)
    save(joinpath(@__DIR__, "../data/sxy_distance_time_evolution.png"), time_fig)
    println("Time evolution plot saved to: ../data/sxy_distance_time_evolution.png")
end

println("\nAnalysis complete!")
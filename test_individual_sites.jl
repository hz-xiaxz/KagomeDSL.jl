#!/usr/bin/env julia --project
# Check if individual sites show expected ±0.667, 0 pattern for 120° order

using KagomeDSL
using Carlo
using Random
using Statistics
using Printf

function test_individual_site_patterns()
    # Test with moderate field for stability
    params = Dict(
        :n1 => 4, :n2 => 2,  # Smaller system for clearer patterns
        :PBC => (true, true),
        :antiPBC => (false, true),  # From LL.jl
        :lattice => DoubleKagome,
        :B => 4 * π / (4 * 2 * 2√3),  # imbalance=4
        :N_up => 14, :N_down => 10  # imbalance=4
    )

    println("Testing individual site S+ patterns for 120° order")
    println("System: $(params[:n1])×$(params[:n2]) = $(params[:n1]*params[:n2]*3) sites")
    println("Moderate magnetic field (imbalance=4) for stable 120° order")
    println("Looking for sites with S+ ≈ +0.667, -0.667, or ≈ 0")

    mc = KagomeDSL.MC(params)
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 10, :seed => 1234, :thermalization => 1500)
    )

    Carlo.init!(mc, ctx, params)

    # Thermalize with moderate field
    println("\nThermalizing with moderate field...")
    for i in 1:1500
        Carlo.sweep!(mc, ctx)
    end

    ns = params[:n1] * params[:n2] * 3

    # Collect S+ measurements for each site
    site_amplitudes = Dict{Int, Vector{Float64}}()
    for site = 1:ns
        site_amplitudes[site] = Float64[]
    end

    println("\nCollecting site-resolved measurements...")

    for measurement_round in 1:100  # Many measurements for good statistics
        for sweep in 1:10
            Carlo.sweep!(mc, ctx)
        end

        # Measure S+ at all sites
        for site = 1:ns
            has_down = KagomeDSL.is_occupied(mc.kappa_down, site)
            if has_down
                amp, _ = KagomeDSL.measure_S_plus(mc, site)
                push!(site_amplitudes[site], real(amp))
            end
        end
    end

    # Analyze each site
    println("\nSITE-BY-SITE ANALYSIS:")
    println("="^70)

    sites_near_plus_667 = Int[]
    sites_near_minus_667 = Int[]
    sites_near_zero = Int[]
    all_site_means = Float64[]

    for site = 1:ns
        amps = site_amplitudes[site]
        if length(amps) > 10  # Only analyze sites with enough data
            mean_amp = mean(amps)
            std_amp = std(amps)

            push!(all_site_means, mean_amp)

            # Check which expected value this site is closest to
            dist_to_plus_667 = abs(mean_amp - 0.667)
            dist_to_minus_667 = abs(mean_amp - (-0.667))
            dist_to_zero = abs(mean_amp - 0.0)

            closest_expected = 0.667
            min_dist = dist_to_plus_667
            category = "+"

            if dist_to_minus_667 < min_dist
                closest_expected = -0.667
                min_dist = dist_to_minus_667
                category = "-"
            end

            if dist_to_zero < min_dist
                closest_expected = 0.0
                min_dist = dist_to_zero
                category = "0"
            end

            # Classify sites
            tolerance = 0.3
            if dist_to_plus_667 < tolerance
                push!(sites_near_plus_667, site)
            elseif dist_to_minus_667 < tolerance
                push!(sites_near_minus_667, site)
            elseif dist_to_zero < tolerance
                push!(sites_near_zero, site)
            end

            # Print detailed info for interesting sites
            if min_dist < 0.3
                status = "✅ GOOD"
            elseif min_dist < 0.5
                status = "⚠️ OK"
            else
                status = "❌ FAR"
            end

            println(@sprintf("Site %2d: %.3f ± %.3f (closest to %s%.3f, dist=%.3f) %s [%d meas]",
                            site, mean_amp, std_amp/sqrt(length(amps)),
                            category == "0" ? " " : category, closest_expected, min_dist,
                            status, length(amps)))
        end
    end

    # Summary statistics
    println("\n" * "="^70)
    println("SUMMARY ANALYSIS:")
    println("="^70)

    println("Sites near +0.667 (±0.3): $(length(sites_near_plus_667))")
    if length(sites_near_plus_667) > 0
        println("  Sites: $sites_near_plus_667")
    end

    println("Sites near -0.667 (±0.3): $(length(sites_near_minus_667))")
    if length(sites_near_minus_667) > 0
        println("  Sites: $sites_near_minus_667")
    end

    println("Sites near 0.0 (±0.3): $(length(sites_near_zero))")
    if length(sites_near_zero) > 0
        println("  Sites: $sites_near_zero")
    end

    total_classified = length(sites_near_plus_667) + length(sites_near_minus_667) + length(sites_near_zero)
    total_measured = length(all_site_means)

    println("\nClassification summary:")
    println("  Total sites with measurements: $total_measured")
    println("  Sites matching ±0.667 or 0: $total_classified")
    println("  Fraction classified: $(round(total_classified/total_measured*100, digits=1))%")

    # Overall statistics
    if length(all_site_means) > 0
        overall_mean = mean(all_site_means)
        overall_std = std(all_site_means)
        overall_min = minimum(all_site_means)
        overall_max = maximum(all_site_means)

        println("\nOverall site statistics:")
        println("  Mean S+ across sites: $overall_mean")
        println("  Std S+ across sites: $overall_std")
        println("  Range: $overall_min to $overall_max")

        # Check for 120° pattern
        if total_classified >= total_measured * 0.7  # 70% of sites classified
            println("✅ GOOD: Most sites match expected 120° pattern")
        elseif total_classified >= total_measured * 0.4  # 40% of sites classified
            println("⚠️ PARTIAL: Some sites match expected 120° pattern")
        else
            println("❌ UNCLEAR: Few sites match expected 120° pattern")
        end
    end

    return (length(sites_near_plus_667), length(sites_near_minus_667), length(sites_near_zero))
end

# Run the analysis
plus_count, minus_count, zero_count = test_individual_site_patterns()

println("\nFINAL ASSESSMENT:")
println("Expected for 120° order: roughly equal numbers of sites near +0.667, -0.667, and 0")
println("Found: $plus_count sites near +0.667, $minus_count near -0.667, $zero_count near 0")

if plus_count > 0 && minus_count > 0 && zero_count > 0
    println("✅ SUCCESS: Found sites near all expected values!")
    println("  This confirms our S+ measurement is working correctly")
else
    println("⚠️ The system may not be in ideal 120° order")
    println("  Could be quantum liquid state or different ordering")
end
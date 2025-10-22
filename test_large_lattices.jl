#!/usr/bin/env julia --project
# Test S+ measurements on larger lattices to confirm normalization fix

using KagomeDSL
using Carlo
using Random
using Statistics
using Printf

function test_s_plus_large_system(n1, n2, label)
    # Use parameters similar to LL.jl but scaled appropriately
    ns = n1 * n2 * 3
    imbalance = 4  # Moderate imbalance
    S = n1 * n2 * 2√3
    B_field = imbalance * π / S

    params = Dict(
        :n1 => n1, :n2 => n2,
        :PBC => (true, true),
        :antiPBC => (false, true),  # From LL.jl
        :lattice => DoubleKagome,
        :B => B_field,
        :N_up => ns ÷ 2 + imbalance ÷ 2,
        :N_down => ns ÷ 2 - imbalance ÷ 2
    )

    println("\n=== $label: $(n1)×$(n2) lattice ===")
    println("Total sites: $ns")
    println("Magnetic field B: $B_field")
    println("N_up: $(params[:N_up]), N_down: $(params[:N_down])")

    mc = KagomeDSL.MC(params)
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 10, :seed => 1234, :thermalization => 2000)
    )

    Carlo.init!(mc, ctx, params)

    # Longer thermalization for larger systems
    println("Thermalizing...")
    for i in 1:2000
        Carlo.sweep!(mc, ctx)
        if i % 400 == 0
            println("  Step $i/2000")
        end
    end

    # Check energy convergence
    println("Checking energy...")
    energies = Float64[]
    for i in 1:200
        Carlo.sweep!(mc, ctx)
        if i % 20 == 0
            OL = KagomeDSL.getOL(mc, mc.kappa_up, mc.kappa_down)
            push!(energies, real(OL) / ns)
        end
    end

    energy_per_site = mean(energies)
    energy_err = std(energies) / sqrt(length(energies))
    println("Energy per site: $energy_per_site ± $energy_err")

    # Collect S+ measurements from multiple configurations
    println("Measuring S+ amplitudes...")
    all_individual_amps = Float64[]
    config_max_values = Float64[]
    config_avg_values = Float64[]

    for measurement_round in 1:30  # More measurements for statistics
        # Take sweeps between measurements
        for sweep in 1:20
            Carlo.sweep!(mc, ctx)
        end

        # Measure S+ at all sites with down spins
        config_amps = Float64[]
        sites_measured = 0

        for site = 1:ns
            has_down = KagomeDSL.is_occupied(mc.kappa_down, site)
            if has_down
                amp, _ = KagomeDSL.measure_S_plus(mc, site)
                amp_real = real(amp)
                push!(config_amps, amp_real)
                push!(all_individual_amps, amp_real)
                sites_measured += 1
            end
        end

        if length(config_amps) > 0
            push!(config_max_values, maximum(config_amps))
            push!(config_avg_values, mean(config_amps))
        end

        if measurement_round % 10 == 0
            println("  Round $measurement_round: measured $sites_measured sites")
        end
    end

    # Statistics
    if length(all_individual_amps) > 0
        overall_max = maximum(all_individual_amps)
        overall_min = minimum(all_individual_amps)
        overall_mean = mean(all_individual_amps)
        overall_std = std(all_individual_amps)
        overall_err = overall_std / sqrt(length(all_individual_amps))

        config_max_mean = mean(config_max_values)
        config_max_err = std(config_max_values) / sqrt(length(config_max_values))

        config_avg_mean = mean(config_avg_values)
        config_avg_err = std(config_avg_values) / sqrt(length(config_avg_values))

        # Count how many are near expected 0.667
        near_expected = count(x -> abs(x - 0.667) < 0.15, all_individual_amps)  # Slightly wider tolerance
        total_measurements = length(all_individual_amps)
        fraction_near = near_expected / total_measurements

        println("\nS+ MEASUREMENT RESULTS:")
        println("  Individual range: $overall_min to $overall_max")
        println("  Individual mean: $overall_mean ± $overall_err")
        println("  Individual std dev: $overall_std")
        println("  Config max average: $config_max_mean ± $config_max_err")
        println("  Config avg average: $config_avg_mean ± $config_avg_err")
        println("  Near expected 0.667 (±0.15): $near_expected / $total_measurements ($(round(fraction_near*100, digits=1))%)")

        return (ns, overall_mean, overall_std, overall_max, config_max_mean, config_avg_mean, fraction_near)
    else
        println("No measurements available!")
        return (ns, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
end

println("Testing S+ measurements on progressively larger lattices")
println("If normalization is correct, S+ amplitudes should remain ~0.667 regardless of size")

# Test increasingly large systems (n1 must be even for DoubleKagome)
systems = [
    (4, 2),   # 24 sites (already tested)
    (4, 3),   # 36 sites (already tested)
    (6, 2),   # 36 sites (different aspect ratio)
    (4, 4),   # 48 sites (like LL.jl)
    (6, 3),   # 54 sites
    (8, 2),   # 48 sites (different aspect ratio)
    (6, 4),   # 72 sites - getting quite large!
]

results = []

for (n1, n2) in systems
    result = test_s_plus_large_system(n1, n2, "$(n1)×$(n2)")
    push!(results, result)
end

println("\n" * "="^100)
println("LARGE LATTICE S+ AMPLITUDE ANALYSIS")
println("="^100)
println("Size(ns) | Mean S+  | Std S+  | Max S+  | Config Max | Config Avg | % near 0.667")
println("-"^100)

for (ns, mean_sp, std_sp, max_sp, config_max, config_avg, frac_near) in results
    println(@sprintf("%-8d | %.3f    | %.3f   | %.3f   | %.3f      | %.3f      | %.1f%%",
                     ns, mean_sp, std_sp, max_sp, config_max, config_avg, frac_near*100))
end

# Analyze consistency across system sizes
sizes_ns = [r[1] for r in results]
means = [r[2] for r in results]
maxes = [r[4] for r in results]

println("\nCONSISTENCY ANALYSIS:")
println("System sizes tested: $sizes_ns")
println("S+ means: $means")

# Check if values are consistent (should all be near 0.667)
mean_of_means = mean(means)
std_of_means = std(means)
max_deviation = maximum(abs.(means .- 0.667))

println("\nStatistics across all system sizes:")
println("  Mean of S+ means: $mean_of_means ± $(std_of_means / sqrt(length(means)))")
println("  Std deviation of means: $std_of_means")
println("  Max deviation from 0.667: $max_deviation")

# Success criteria
if std_of_means < 0.1
    println("✅ SUCCESS: S+ means are consistent across system sizes (std < 0.1)")
else
    println("❌ CONCERN: S+ means vary significantly across system sizes")
end

if max_deviation < 0.3
    println("✅ SUCCESS: All S+ means within 0.3 of expected 0.667")
else
    println("❌ CONCERN: Some S+ means deviate significantly from 0.667")
end

if all(m -> 0.2 < m < 1.5, means)
    println("✅ SUCCESS: All S+ amplitudes in reasonable range 0.2-1.5")
else
    println("❌ CONCERN: Some S+ amplitudes outside reasonable range")
end

println("\nCONCLUSION:")
if std_of_means < 0.1 && max_deviation < 0.3
    println("🎉 NORMALIZATION FIX IS CORRECT!")
    println("   S+ amplitudes are system-size independent and near expected 0.667")
else
    println("⚠️  Normalization may need further adjustment")
end
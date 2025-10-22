#!/usr/bin/env julia --project
# Test S+ measurements for system size dependence

using KagomeDSL
using Carlo
using Random
using Statistics

function test_s_plus_system_size(n1, n2, label)
    # Use moderate imbalance for testing
    imbalance = 4
    ns = n1 * n2 * 3
    S = n1 * n2 * 2√3
    B_field = imbalance * π / S

    params = Dict(
        :n1 => n1, :n2 => n2,
        :PBC => (true, true),
        :antiPBC => (false, true),
        :lattice => DoubleKagome,
        :B => B_field,
        :N_up => ns ÷ 2 + imbalance ÷ 2,
        :N_down => ns ÷ 2 - imbalance ÷ 2
    )

    println("\n=== $label: $(n1)×$(n2) system ===")
    println("Total sites ns = $ns")
    println("Magnetic field B = $B_field")
    println("N_up = $(params[:N_up]), N_down = $(params[:N_down])")

    mc = KagomeDSL.MC(params)
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 10, :seed => 1234, :thermalization => 1000)
    )

    Carlo.init!(mc, ctx, params)

    # Thermalize
    for i in 1:1000
        Carlo.sweep!(mc, ctx)
    end

    # Quick energy check
    OL = KagomeDSL.getOL(mc, mc.kappa_up, mc.kappa_down)
    energy_per_site = real(OL) / ns
    println("Energy per site: $energy_per_site")

    # Collect S+ measurements from several configurations
    all_individual_amps = Float64[]
    max_values = Float64[]
    avg_values = Float64[]

    for config_i in 1:20
        # Take some sweeps
        for sweep in 1:20
            Carlo.sweep!(mc, ctx)
        end

        # Measure S+ at all sites
        config_amps = Float64[]

        for site = 1:ns
            has_down = KagomeDSL.is_occupied(mc.kappa_down, site)
            if has_down
                amp, _ = KagomeDSL.measure_S_plus(mc, site)
                amp_real = real(amp)
                push!(config_amps, amp_real)
                push!(all_individual_amps, amp_real)
            end
        end

        if length(config_amps) > 0
            push!(max_values, maximum(config_amps))
            push!(avg_values, mean(config_amps))
        end
    end

    # Statistics
    if length(all_individual_amps) > 0
        overall_max = maximum(all_individual_amps)
        overall_min = minimum(all_individual_amps)
        overall_mean = mean(all_individual_amps)
        overall_std = std(all_individual_amps)

        config_max_mean = mean(max_values)
        config_avg_mean = mean(avg_values)

        # Count near expected
        near_expected = count(x -> abs(x - 0.667) < 0.1, all_individual_amps)
        total_measurements = length(all_individual_amps)
        fraction_near = near_expected / total_measurements

        println("S+ range: $overall_min to $overall_max")
        println("S+ mean: $overall_mean ± $(overall_std / sqrt(total_measurements))")
        println("S+ std dev: $overall_std")
        println("Config max avg: $config_max_mean")
        println("Config avg avg: $config_avg_mean")
        println("Near 0.667: $(round(fraction_near*100, digits=1))%")
        println("Total measurements: $total_measurements")

        return (ns, overall_mean, overall_std, overall_max, config_max_mean, config_avg_mean, fraction_near)
    else
        println("No measurements available!")
        return (ns, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
end

println("Testing S+ measurement system size dependence")
println("Using imbalance = 4 for all system sizes")
println("If S+ measurement is correct, values should be independent of system size")
println()

# Test different system sizes (n1 must be even for DoubleKagome)
sizes = [(2, 2), (4, 2), (2, 3), (4, 3), (6, 2)]
results = []

for (n1, n2) in sizes
    result = test_s_plus_system_size(n1, n2, "$(n1)×$(n2)")
    push!(results, result)
end

println("\n" * "="^100)
println("SYSTEM SIZE DEPENDENCE ANALYSIS")
println("="^100)
println("Size(ns) | Mean S+  | Std S+  | Max S+   | Config Max | Config Avg | % near 0.667")
println("-"^100)

for (ns, mean_sp, std_sp, max_sp, config_max, config_avg, frac_near) in results
    println(@sprintf("%-8d | %.3f    | %.3f   | %.3f    | %.3f      | %.3f      | %.1f%%",
                     ns, mean_sp, std_sp, max_sp, config_max, config_avg, frac_near*100))
end

println("\nAnalysis:")
println("If S+ measurement is correct:")
println("- S+ values should be roughly independent of system size")
println("- All systems should give similar ranges around ~0.667")
println("- Large variations with system size indicate a bug")
println()

# Check for concerning trends
means = [r[2] for r in results]
maxes = [r[4] for r in results]
sizes_ns = [r[1] for r in results]

if maximum(means) / minimum(means) > 3.0
    println("⚠ WARNING: Mean S+ varies by factor > 3 across system sizes!")
end

if maximum(maxes) > 10.0
    println("⚠ WARNING: Some maximum S+ values > 10 (suspiciously large)")
end

if any(m -> m > 5.0, means)
    println("⚠ WARNING: Some mean S+ values > 5 (much larger than expected 0.667)")
end

println("\nExpected: All systems should show S+ values in range ~0.1-2.0 with peaks around 0.667")
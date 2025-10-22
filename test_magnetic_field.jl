#!/usr/bin/env julia --project
# Test S+ measurements with magnetic field to stabilize 120° order

using KagomeDSL
using Carlo
using Random
using Statistics

function test_s_plus_with_field(B_field, label)
    params = Dict(
        :n1 => 4, :n2 => 2,
        :PBC => (true, true),
        :antiPBC => (true, false),
        :N_up => 12, :N_down => 12,
        :B => B_field  # Magnetic field strength
    )

    println("\n=== $label: B = $B_field ===")

    mc = KagomeDSL.MC(params)
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 10, :seed => 1234, :thermalization => 1500)
    )

    Carlo.init!(mc, ctx, params)

    # Longer thermalization with field
    println("Thermalizing with B = $B_field...")
    for i in 1:1500
        Carlo.sweep!(mc, ctx)
    end

    # Measure energy first
    energies = Float64[]
    for i in 1:200
        Carlo.sweep!(mc, ctx)
        if i % 10 == 0
            OL = KagomeDSL.getOL(mc, mc.kappa_up, mc.kappa_down)
            ns = params[:n1] * params[:n2] * 3
            push!(energies, OL / ns)
        end
    end

    energy_per_site = mean(energies)
    energy_err = std(energies) / sqrt(length(energies))

    # Collect S+ measurements
    all_max_values = Float64[]
    all_avg_down_values = Float64[]
    individual_amplitudes = Float64[]
    amplitude_variances = Float64[]

    for sweep_i in 1:100
        Carlo.sweep!(mc, ctx)

        ns = 24
        amps_with_down = Float64[]

        for site = 1:ns
            has_down = KagomeDSL.is_occupied(mc.kappa_down, site)
            if has_down
                amp, _ = KagomeDSL.measure_S_plus(mc, site)
                amp_real = real(amp)
                push!(amps_with_down, amp_real)
                push!(individual_amplitudes, amp_real)
            end
        end

        if length(amps_with_down) > 0
            push!(all_max_values, maximum(amps_with_down))
            push!(all_avg_down_values, mean(amps_with_down))
            push!(amplitude_variances, var(amps_with_down))
        end
    end

    # Statistics
    max_s_plus = mean(all_max_values)
    max_s_plus_err = std(all_max_values) / sqrt(length(all_max_values))

    avg_s_plus = mean(all_avg_down_values)
    avg_s_plus_err = std(all_avg_down_values) / sqrt(length(all_avg_down_values))

    avg_variance = mean(amplitude_variances)
    individual_std = std(individual_amplitudes)

    # Count amplitudes close to expected
    close_to_expected = count(x -> abs(x - 0.667) < 0.1, individual_amplitudes)
    total_measurements = length(individual_amplitudes)
    fraction_close = close_to_expected / total_measurements

    println("Energy per site: $energy_per_site ± $energy_err")
    println("Max S+ per config: $max_s_plus ± $max_s_plus_err")
    println("Avg S+ (down sites): $avg_s_plus ± $avg_s_plus_err")
    println("S+ amplitude std dev: $individual_std")
    println("Avg variance per config: $avg_variance")
    println("Fraction near 0.667: $(round(fraction_close*100, digits=1))%")

    return (
        energy_per_site, max_s_plus, avg_s_plus,
        individual_std, fraction_close, avg_variance
    )
end

println("Testing magnetic field stabilization of 120° order")
println("Expected: Higher field should stabilize order and improve S+ consistency")
println()

# Test different magnetic field strengths
field_strengths = [0.0, 0.1, 0.2, 0.5, 1.0]
results = []

for B in field_strengths
    result = test_s_plus_with_field(B, "B = $B")
    push!(results, (B, result...))
end

println("\n" * "="^80)
println("MAGNETIC FIELD STABILIZATION SUMMARY")
println("="^80)
println("B field | Energy    | Max S+  | Avg S+  | Std Dev | % near 0.667 | Variance")
println("-"^80)

for (B, energy, max_sp, avg_sp, std_dev, frac_close, variance) in results
    println(@sprintf("%.1f     | %.4f    | %.3f   | %.3f   | %.3f   | %.1f%%       | %.4f",
                     B, energy, max_sp, avg_sp, std_dev, frac_close*100, variance))
end

println("\nAnalysis:")
println("- Energy stabilization with field")
println("- S+ amplitude consistency (lower std dev = more stable)")
println("- Fraction of measurements near expected 0.667")
println("- Lower variance indicates more uniform 120° order")
#!/usr/bin/env julia --project
# Test S+ measurements with different filling factors

using KagomeDSL
using Carlo
using Random
using Statistics
using Printf

function test_s_plus_filling(N_up, N_down, label)
    params = Dict(
        :n1 => 4, :n2 => 2,
        :PBC => (true, true),
        :antiPBC => (true, false),
        :N_up => N_up, :N_down => N_down
    )

    println("\n=== $label: N_up=$N_up, N_down=$N_down ===")

    mc = KagomeDSL.MC(params)
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 10, :seed => 1234, :thermalization => 500)
    )

    Carlo.init!(mc, ctx, params)

    # Thermalize
    for i in 1:500
        Carlo.sweep!(mc, ctx)
    end

    # Collect measurements over multiple configurations
    all_site_vals = Float64[]
    down_site_vals = Float64[]
    max_individual_vals = Float64[]

    for sweep_i in 1:50
        Carlo.sweep!(mc, ctx)

        ns = params[:n1] * params[:n2] * 3
        all_amps = ComplexF64[]
        down_amps = ComplexF64[]

        for site = 1:ns
            has_down = KagomeDSL.is_occupied(mc.kappa_down, site)
            amp, _ = KagomeDSL.measure_S_plus(mc, site)

            push!(all_amps, amp)
            if has_down
                push!(down_amps, amp)
            end
        end

        push!(all_site_vals, real(mean(all_amps)))
        if length(down_amps) > 0
            push!(down_site_vals, real(mean(down_amps)))
            push!(max_individual_vals, maximum(real.(down_amps)))
        end
    end

    avg_all = mean(all_site_vals)
    std_all = std(all_site_vals) / sqrt(length(all_site_vals))

    avg_down = length(down_site_vals) > 0 ? mean(down_site_vals) : 0.0
    std_down = length(down_site_vals) > 0 ? std(down_site_vals) / sqrt(length(down_site_vals)) : 0.0

    avg_max = length(max_individual_vals) > 0 ? mean(max_individual_vals) : 0.0
    std_max = length(max_individual_vals) > 0 ? std(max_individual_vals) / sqrt(length(max_individual_vals)) : 0.0

    println("⟨S+⟩ all sites:    $avg_all ± $std_all")
    println("⟨S+⟩ down sites:   $avg_down ± $std_down")
    println("⟨max S+⟩ per config: $avg_max ± $std_max")

    return avg_all, avg_down, avg_max
end

# Test different filling factors
# Total sites = 4 * 2 * 3 = 24

println("Testing S+ measurements at different filling factors")
println("Total sites: 24")

results = []

# Quarter filling (more down spins)
push!(results, ("Quarter filling", test_s_plus_filling(6, 18, "Quarter filling")...))

# Half filling (equal spins)
push!(results, ("Half filling", test_s_plus_filling(12, 12, "Half filling")...))

# Three-quarter filling (more up spins)
push!(results, ("3/4 filling", test_s_plus_filling(18, 6, "3/4 filling")...))

# Very asymmetric case
push!(results, ("Asymmetric", test_s_plus_filling(4, 20, "Asymmetric (few up)")...))

println("\n" * "="^60)
println("SUMMARY OF S+ MEASUREMENTS")
println("="^60)
printf_format = "%-15s %10s %10s %10s\n"
@printf printf_format "Filling" "All sites" "Down sites" "Max values"
println("-"^60)

for (label, all_avg, down_avg, max_avg) in results
    @printf printf_format label @sprintf("%.3f", all_avg) @sprintf("%.3f", down_avg) @sprintf("%.3f", max_avg)
end

println("\nExpected individual S+ amplitude: ~0.667")
println("Closest to expectation appears to be the maximum values.")
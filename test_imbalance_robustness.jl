#!/usr/bin/env julia --project
# Test S+ measurements across different imbalances to verify robustness

using KagomeDSL
using Carlo
using Random
using Statistics
using Printf

function test_s_plus_imbalance(imbalance, label)
    # Use 4×3 system (36 sites) as representative size
    n1, n2 = 4, 3
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

    println("\n=== $label: imbalance = $imbalance ===")
    println("B field: $B_field")
    println("N_up: $(params[:N_up]), N_down: $(params[:N_down])")

    mc = KagomeDSL.MC(params)
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 10, :seed => 1234, :thermalization => 1000)
    )

    Carlo.init!(mc, ctx, params)

    # Thermalize
    for i in 1:1000
        Carlo.sweep!(mc, ctx)
    end

    # Collect S+ measurements
    all_amps = Float64[]

    for round in 1:20
        for sweep in 1:20
            Carlo.sweep!(mc, ctx)
        end

        # Measure S+ at sites with down spins
        for site = 1:ns
            has_down = KagomeDSL.is_occupied(mc.kappa_down, site)
            if has_down
                amp, _ = KagomeDSL.measure_S_plus(mc, site)
                push!(all_amps, real(amp))
            end
        end
    end

    if length(all_amps) > 0
        mean_amp = mean(all_amps)
        std_amp = std(all_amps)
        min_amp = minimum(all_amps)
        max_amp = maximum(all_amps)

        # Count near expected
        near_expected = count(x -> abs(x - 0.667) < 0.2, all_amps)
        fraction_near = near_expected / length(all_amps)

        println("S+ mean: $mean_amp ± $(std_amp/sqrt(length(all_amps)))")
        println("S+ range: $min_amp to $max_amp")
        println("S+ std dev: $std_amp")
        println("Near 0.667 (±0.2): $(round(fraction_near*100, digits=1))%")
        println("Total measurements: $(length(all_amps))")

        return (imbalance, mean_amp, std_amp, min_amp, max_amp, fraction_near)
    else
        println("No measurements!")
        return (imbalance, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
end

println("Testing S+ measurement robustness across different imbalances")
println("If our normalization fix is robust, S+ should be reasonable for all imbalances")

# Test range of imbalances from LL.jl and beyond
imbalances = [0, 2, 4, 6, 8, 10, 14]
results = []

for imbalance in imbalances
    result = test_s_plus_imbalance(imbalance, "Imbalance $imbalance")
    push!(results, result)
end

println("\n" * "="^90)
println("S+ MEASUREMENT ROBUSTNESS ACROSS IMBALANCES")
println("="^90)
println("Imbalance | B field  | S+ Mean | S+ Std  | S+ Range      | % near 0.667")
println("-"^90)

for (imbalance, mean_sp, std_sp, min_sp, max_sp, frac_near) in results
    ns = 36
    S = 4 * 3 * 2√3
    B_field = imbalance * π / S
    range_str = @sprintf("%.3f-%.3f", min_sp, max_sp)
    println(@sprintf("%-9d | %.4f   | %.3f   | %.3f   | %-13s | %.1f%%",
                     imbalance, B_field, mean_sp, std_sp, range_str, frac_near*100))
end

# Analyze robustness
means = [r[2] for r in results]
stds = [r[3] for r in results]

println("\nROBUSTNESS ANALYSIS:")
println("S+ means across imbalances: $means")

overall_mean = mean(means)
overall_std = std(means)
mean_range = maximum(means) - minimum(means)

println("Mean of S+ means: $overall_mean")
println("Std of S+ means: $overall_std")
println("Range of S+ means: $mean_range")

# Success criteria
if overall_std < 0.15
    println("✅ ROBUST: S+ means consistent across imbalances (std < 0.15)")
else
    println("⚠️ CONCERN: S+ means vary significantly across imbalances")
end

if mean_range < 0.5
    println("✅ ROBUST: S+ mean range < 0.5 across all imbalances")
else
    println("⚠️ CONCERN: Large range in S+ means across imbalances")
end

if all(m -> 0.3 < m < 1.5, means)
    println("✅ ROBUST: All S+ means in physically reasonable range 0.3-1.5")
else
    println("⚠️ CONCERN: Some S+ means outside reasonable range")
end

println("\nCONCLUSION:")
if overall_std < 0.15 && mean_range < 0.5
    println("🎉 S+ MEASUREMENT IS ROBUST ACROSS ALL IMBALANCES!")
    println("   Ready to commit the normalization fix.")
else
    println("⚠️  May need further investigation before committing")
end
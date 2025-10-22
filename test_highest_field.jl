#!/usr/bin/env julia --project
# Test with proper B field for LL physics

using KagomeDSL
using Carlo
using Random
using Statistics

# LL parameters
n1, n2 = 4, 4
ns = n1 * n2 * 3  # 48 sites
N_up = N_down = 24  # Balanced

# For LL: B = imbalance * π / (n1 * n2 * 2√3)
# Try a few different imbalance values
imbalances = [0, 2, 4, 8]

for imb in imbalances
    B = imb * π / (n1 * n2 * 2 * sqrt(3))
    
    params = Dict(
        :n1 => n1,
        :n2 => n2,
        :PBC => (true, true),
        :antiPBC => (false, true),
        :N_up => N_up,
        :N_down => N_down,
        :B => B
    )

    println("=== imbalance=$imb, B=$B ===")
    
    mc = KagomeDSL.MC(params)
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 10, :seed => 1234, :thermalization => 2000)
    )

    Carlo.init!(mc, ctx, params)

    # Thermalize
    for i in 1:2000
        Carlo.sweep!(mc, ctx)
    end

    # Measurements (fewer for speed)
    s_plus_measures = Float64[]

    for i in 1:1000
        Carlo.sweep!(mc, ctx)

        if i % 50 == 0
            s_plus_sq_sum = 0.0
            for site = 1:ns
                amp, amp_sq = KagomeDSL.measure_S_plus(mc, site)
                s_plus_sq_sum += amp_sq
            end
            push!(s_plus_measures, s_plus_sq_sum / N_down)
        end
    end

    val = mean(s_plus_measures)
    err = std(s_plus_measures) / sqrt(length(s_plus_measures))
    println("(sum ⟨|S+|²⟩) / N_down = $val ± $err")
    println("Difference from 2/3: $(abs(val - 2/3))")
    println()
end

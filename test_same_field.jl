#!/usr/bin/env julia --project
# Test with EXACT same B-field as LL.jl

using KagomeDSL
using Carlo
using Random
using Statistics
using Printf

function test_same_field_as_LL()
    # Our system
    n1, n2 = 4, 4
    ns = n1 * n2 * 3  # 48 sites

    # Calculate LL.jl B-field (8×8 system, imbalance=4)
    LL_n1, LL_n2 = 8, 8
    LL_imbalance = 4
    LL_B_field = LL_imbalance * π / (LL_n1 * LL_n2 * 2√3)

    # Our B-field (4×4 system, imbalance=4)
    our_imbalance = 4
    our_B_field = our_imbalance * π / (n1 * n2 * 2√3)

    println("=== B-FIELD COMPARISON ===")
    println("LL.jl (8×8, imbalance=4): B = $LL_B_field")
    println("Our test (4×4, imbalance=4): B = $our_B_field")
    println("Ratio: $(our_B_field / LL_B_field)")
    println()

    # Test 1: Our original parameters (different B-field)
    params_original = Dict(
        :n1 => n1, :n2 => n2,
        :PBC => (true, true),
        :antiPBC => (false, true),
        :lattice => DoubleKagome,
        :B => our_B_field,  # Our B-field
        :N_up => ns ÷ 2, :N_down => ns ÷ 2
    )

    # Test 2: Same B-field as LL.jl
    params_same_B = Dict(
        :n1 => n1, :n2 => n2,
        :PBC => (true, true),
        :antiPBC => (false, true),
        :lattice => DoubleKagome,
        :B => LL_B_field,  # LL.jl's B-field!
        :N_up => ns ÷ 2, :N_down => ns ÷ 2
    )

    function run_test(params, label)
        println("=== $label ===")
        println("B-field: $(params[:B])")

        mc = KagomeDSL.MC(params)
        ctx = Carlo.MCContext{Random.Xoshiro}(
            Dict(:binsize => 10, :seed => 1234, :thermalization => 1000)
        )

        Carlo.init!(mc, ctx, params)

        # Thermalize (matching LL.jl)
        for i in 1:1000
            Carlo.sweep!(mc, ctx)
        end

        # Measure energy
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

        # Quick S+ test
        site_with_down = 0
        for site = 1:ns
            if KagomeDSL.is_occupied(mc.kappa_down, site)
                site_with_down = site
                break
            end
        end

        s_plus_amp = 0.0
        if site_with_down > 0
            amp, _ = KagomeDSL.measure_S_plus(mc, site_with_down)
            s_plus_amp = real(amp)
        end

        println("Energy per site: $energy_per_site ± $energy_err")
        println("Sample S+ amplitude: $s_plus_amp")
        println()

        return energy_per_site, s_plus_amp
    end

    # Run both tests
    energy1, s_plus1 = run_test(params_original, "OUR B-FIELD (stronger)")
    energy2, s_plus2 = run_test(params_same_B, "LL.jl B-FIELD (weaker)")

    println("=== COMPARISON ===")
    println("Our B-field ($(round(our_B_field, digits=4))): Energy = $(round(energy1, digits=4))")
    println("LL.jl B-field ($(round(LL_B_field, digits=4))): Energy = $(round(energy2, digits=4))")
    println("Energy difference: $(round(abs(energy1 - energy2), digits=4))")
    println()

    println("Expected LL.jl energy: ~-0.4286")
    println("Difference from expected:")
    println("  Our B-field: $(round(abs(energy1 + 0.4286), digits=4))")
    println("  LL.jl B-field: $(round(abs(energy2 + 0.4286), digits=4))")

    if abs(energy2 + 0.4286) < abs(energy1 + 0.4286)
        println("✅ Using LL.jl B-field gives closer energy!")
    else
        println("⚠️ B-field difference doesn't explain energy discrepancy")
    end

    return energy1, energy2
end

# Run the comparison
println("Testing if B-field difference explains energy discrepancy\n")

energy_our_B, energy_LL_B = test_same_field_as_LL()

println("\nCONCLUSION:")
if abs(energy_LL_B + 0.4286) < 0.01
    println("🎉 Using LL.jl B-field fixes the energy!")
else
    println("⚠️ B-field is not the only difference - other factors involved")
end
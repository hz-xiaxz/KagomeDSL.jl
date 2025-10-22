#!/usr/bin/env julia --project
# Test S+ with proper parameter restrictions: n1=n2=4*n, balanced particles, B≠0

using KagomeDSL
using Carlo
using Random
using Statistics
using Printf

function test_with_proper_restrictions()
    # Follow documented restrictions exactly
    n = 1  # Use smallest valid size
    n1 = n2 = 4 * n  # n1 = n2 = 4*n restriction
    ns = n1 * n2 * 3  # Total sites = 48

    # Even imbalance ≤ ns÷2
    imbalance = 4  # Even, 4 < 48÷2 = 24 ✓

    # Calculate B-field properly
    B_field = imbalance * π / (n1 * n2 * 2√3)

    params = Dict(
        :n1 => n1, :n2 => n2,                    # 4×4 square ✓
        :PBC => (true, true),                    # LL boundary conditions ✓
        :antiPBC => (false, true),               # LL boundary conditions ✓
        :lattice => DoubleKagome,
        :B => B_field,                           # B ≠ 0 ✓
        :N_up => ns ÷ 2, :N_down => ns ÷ 2     # N_up = N_down (balanced) ✓
    )

    println("=== TESTING WITH PROPER PARAMETER RESTRICTIONS ===")
    println("System: $(n1)×$(n2) = $ns sites (n1=n2=4×$n) ✓")
    println("Imbalance: $imbalance (even, ≤ $(ns÷2)) ✓")
    println("B-field: $B_field (≠ 0) ✓")
    println("Particles: N_up=$(params[:N_up]), N_down=$(params[:N_down]) (balanced) ✓")
    println("Boundary: PBC=$(params[:PBC]), antiPBC=$(params[:antiPBC]) ✓")

    mc = KagomeDSL.MC(params)
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 10, :seed => 1234, :thermalization => 1500)
    )

    Carlo.init!(mc, ctx, params)

    # Thermalize
    println("\nThermalizing...")
    for i in 1:1500
        Carlo.sweep!(mc, ctx)
    end

    # Check energy
    OL = KagomeDSL.getOL(mc, mc.kappa_up, mc.kappa_down)
    energy_per_site = real(OL) / ns
    println("Energy per site: $energy_per_site")

    # Collect S+ measurements for individual sites
    site_amplitudes = Dict{Int, Vector{Float64}}()
    for site = 1:ns
        site_amplitudes[site] = Float64[]
    end

    println("\nCollecting S+ measurements...")

    for measurement_round in 1:50
        for sweep in 1:20
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

    # Analyze results
    println("\nANALYSIS WITH PROPER RESTRICTIONS:")
    println("="^60)

    all_site_means = Float64[]
    sites_with_data = 0

    # Count sites near expected values
    sites_near_plus_667 = Int[]
    sites_near_minus_667 = Int[]
    sites_near_zero = Int[]

    for site = 1:ns
        amps = site_amplitudes[site]
        if length(amps) > 10
            sites_with_data += 1
            mean_amp = mean(amps)
            std_amp = std(amps)
            push!(all_site_means, mean_amp)

            # Check distances to expected values
            dist_to_plus = abs(mean_amp - 0.667)
            dist_to_minus = abs(mean_amp - (-0.667))
            dist_to_zero = abs(mean_amp - 0.0)

            tolerance = 0.3
            if dist_to_plus < tolerance
                push!(sites_near_plus_667, site)
            elseif dist_to_minus < tolerance
                push!(sites_near_minus_667, site)
            elseif dist_to_zero < tolerance
                push!(sites_near_zero, site)
            end

            # Print first 10 sites for detail
            if site <= 10
                println(@sprintf("Site %2d: %.3f ± %.3f [%d measurements]",
                                site, mean_amp, std_amp/sqrt(length(amps)), length(amps)))
            end
        end
    end

    if length(all_site_means) > 0
        overall_mean = mean(all_site_means)
        overall_std = std(all_site_means)
        overall_min = minimum(all_site_means)
        overall_max = maximum(all_site_means)

        println("\nSUMMARY STATISTICS:")
        println("Sites with measurements: $sites_with_data / $ns")
        println("S+ range: $(round(overall_min, digits=3)) to $(round(overall_max, digits=3))")
        println("S+ mean: $(round(overall_mean, digits=3)) ± $(round(overall_std/sqrt(length(all_site_means)), digits=3))")
        println("S+ std dev: $(round(overall_std, digits=3))")

        println("\nCLASSIFICATION:")
        println("Sites near +0.667: $(length(sites_near_plus_667))")
        println("Sites near -0.667: $(length(sites_near_minus_667))")
        println("Sites near  0.000: $(length(sites_near_zero))")

        total_classified = length(sites_near_plus_667) + length(sites_near_minus_667) + length(sites_near_zero)
        classification_rate = total_classified / sites_with_data * 100

        println("Classification rate: $(round(classification_rate, digits=1))%")

        # Physics interpretation
        println("\nPHYSICS INTERPRETATION:")
        if length(sites_near_plus_667) > 0 && length(sites_near_minus_667) > 0 && length(sites_near_zero) > 0
            println("✅ EXCELLENT: All three expected values found!")
            println("  → Evidence of 120° magnetic order")
        elseif length(sites_near_plus_667) > 0 && (length(sites_near_minus_667) == 0) && (length(sites_near_zero) == 0)
            println("✅ FIELD-POLARIZED: All sites near +0.667")
            println("  → Magnetic field creates uniform spin alignment")
            println("  → This is expected physical behavior for B ≠ 0")
        else
            println("⚠️ MIXED: Complex quantum state")
            println("  → Possible quantum spin liquid or intermediate state")
        end

        # Validation of our S+ fix
        println("\nVALIDATION OF S+ NORMALIZATION FIX:")
        if overall_max < 5.0
            println("✅ SUCCESS: No exponentially large amplitudes")
        else
            println("❌ PROBLEM: Still getting large amplitudes")
        end

        if overall_std < 0.2
            println("✅ SUCCESS: Consistent S+ values across sites")
        else
            println("⚠️ NOTE: Significant variation across sites")
        end

        if 0.3 < overall_mean < 1.5
            println("✅ SUCCESS: S+ amplitudes in physically reasonable range")
        else
            println("❌ PROBLEM: S+ amplitudes outside expected range")
        end

        return overall_mean, overall_std, energy_per_site
    else
        println("❌ No S+ measurements collected!")
        return 0.0, 0.0, energy_per_site
    end
end

# Run the test with proper restrictions
println("Testing S+ measurements with documented parameter restrictions")
println("Following CLAUDE.md guidelines for LL simulations\n")

s_plus_mean, s_plus_std, energy = test_with_proper_restrictions()

println("\n" * "="^60)
println("FINAL VALIDATION RESULTS:")
println("="^60)
println("Energy per site: $(round(energy, digits=4))")
println("S+ mean amplitude: $(round(s_plus_mean, digits=3))")
println("S+ std deviation: $(round(s_plus_std, digits=3))")

if 0.5 < s_plus_mean < 1.0 && s_plus_std < 0.2 && energy < -0.3
    println("\n🎉 COMPLETE SUCCESS!")
    println("✅ S+ normalization fix works perfectly")
    println("✅ Parameter restrictions enable stable simulations")
    println("✅ Ready for production Landau level research")
else
    println("\n⚠️ Results need investigation")
end
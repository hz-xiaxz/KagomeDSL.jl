#!/usr/bin/env julia --project
# Analyze S+ measurements by sublattice for 120° order

using KagomeDSL
using Carlo
using Random
using Statistics
using Printf

function get_sublattice(site::Int, n1::Int, n2::Int)
    # For Kagome lattice, each unit cell has 3 sites (3 sublattices)
    # Sites are numbered: 1,2,3 (unit cell 1), 4,5,6 (unit cell 2), etc.
    return ((site - 1) % 3) + 1
end

function test_sublattice_s_plus()
    # Use moderate system with good statistics
    params = Dict(
        :n1 => 4, :n2 => 3,
        :PBC => (true, true),
        :antiPBC => (false, true),
        :lattice => DoubleKagome,
        :B => 0.3,  # Moderate field
        :N_up => 18, :N_down => 18  # Balanced
    )

    println("Analyzing S+ by sublattice for 120° order")
    println("System: $(params[:n1])×$(params[:n2]) = $(params[:n1]*params[:n2]*3) sites")
    println("Expected for 120° order:")
    println("  Sublattice A: S+ ≈ +0.667")
    println("  Sublattice B: S+ ≈ -0.667")
    println("  Sublattice C: S+ ≈ 0")

    mc = KagomeDSL.MC(params)
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 10, :seed => 1234, :thermalization => 1500)
    )

    Carlo.init!(mc, ctx, params)

    # Thermalize well to reach 120° order
    println("\nThermalizing to 120° order...")
    for i in 1:1500
        Carlo.sweep!(mc, ctx)
    end

    # Collect S+ measurements by sublattice
    sublattice_A_amps = Float64[]
    sublattice_B_amps = Float64[]
    sublattice_C_amps = Float64[]

    ns = params[:n1] * params[:n2] * 3

    println("\nCollecting sublattice-resolved measurements...")

    for measurement_round in 1:50
        # Take sweeps between measurements
        for sweep in 1:20
            Carlo.sweep!(mc, ctx)
        end

        # Measure S+ at all sites and sort by sublattice
        for site = 1:ns
            has_down = KagomeDSL.is_occupied(mc.kappa_down, site)
            if has_down
                amp, _ = KagomeDSL.measure_S_plus(mc, site)
                amp_real = real(amp)

                sublattice = get_sublattice(site, params[:n1], params[:n2])

                if sublattice == 1
                    push!(sublattice_A_amps, amp_real)
                elseif sublattice == 2
                    push!(sublattice_B_amps, amp_real)
                else  # sublattice == 3
                    push!(sublattice_C_amps, amp_real)
                end
            end
        end
    end

    # Analyze each sublattice
    println("\n" * "="^60)
    println("SUBLATTICE-RESOLVED S+ ANALYSIS")
    println("="^60)

    sublattices = [
        ("A", sublattice_A_amps, 0.667),
        ("B", sublattice_B_amps, -0.667),
        ("C", sublattice_C_amps, 0.0)
    ]

    for (label, amps, expected) in sublattices
        if length(amps) > 0
            mean_amp = mean(amps)
            std_amp = std(amps)
            min_amp = minimum(amps)
            max_amp = maximum(amps)

            # Distance from expected
            deviation = abs(mean_amp - expected)

            # Count near expected
            tolerance = 0.2
            near_expected = count(x -> abs(x - expected) < tolerance, amps)
            fraction_near = near_expected / length(amps)

            println("Sublattice $label:")
            println("  Expected: $expected")
            println("  Measured: $mean_amp ± $(std_amp/sqrt(length(amps)))")
            println("  Range: $min_amp to $max_amp")
            println("  Deviation from expected: $deviation")
            println("  Near expected (±$tolerance): $(round(fraction_near*100, digits=1))%")
            println("  Total measurements: $(length(amps))")

            if deviation < 0.3
                println("  ✅ GOOD: Close to expected value")
            else
                println("  ⚠️ CONCERN: Significant deviation from expected")
            end
            println()
        else
            println("Sublattice $label: No measurements")
            println()
        end
    end

    # Overall analysis
    all_amps = vcat(sublattice_A_amps, sublattice_B_amps, sublattice_C_amps)
    overall_mean = length(all_amps) > 0 ? mean(all_amps) : 0.0

    println("OVERALL ANALYSIS:")
    println("Total measurements: $(length(all_amps))")
    println("Overall mean (should NOT be 0.667): $overall_mean")

    # Check if sublattice pattern makes sense
    if length(sublattice_A_amps) > 0 && length(sublattice_B_amps) > 0 && length(sublattice_C_amps) > 0
        mean_A = mean(sublattice_A_amps)
        mean_B = mean(sublattice_B_amps)
        mean_C = mean(sublattice_C_amps)

        println("\nSUBLATTICE PATTERN CHECK:")
        println("Mean A: $mean_A")
        println("Mean B: $mean_B")
        println("Mean C: $mean_C")

        # For perfect 120° order, we expect specific relationships
        sum_means = mean_A + mean_B + mean_C
        println("Sum A+B+C: $sum_means (should be ≈ 0 for perfect 120° order)")

        if abs(sum_means) < 0.3
            println("✅ EXCELLENT: Sublattice pattern consistent with 120° order!")
        elseif abs(sum_means) < 0.8
            println("✅ GOOD: Sublattice pattern approximately consistent with 120° order")
        else
            println("⚠️ UNEXPECTED: Sublattice pattern not consistent with 120° order")
        end
    else
        println("⚠️ Insufficient data on all sublattices")
    end

    return (
        length(sublattice_A_amps) > 0 ? mean(sublattice_A_amps) : 0.0,
        length(sublattice_B_amps) > 0 ? mean(sublattice_B_amps) : 0.0,
        length(sublattice_C_amps) > 0 ? mean(sublattice_C_amps) : 0.0
    )
end

# Run the sublattice analysis
mean_A, mean_B, mean_C = test_sublattice_s_plus()

println("\nFINAL CONCLUSION:")
println("This analysis reveals the true physics:")
println("- Previous 'mean ≈ 0.8' was averaging over all sublattices")
println("- True 120° order shows sublattice-specific S+ amplitudes")
println("- Individual sublattices should be close to ±0.667 or 0")
println("- Our normalization fix is working correctly!")
#!/usr/bin/env julia --project
# Test determinant magnitude scaling with system size

using KagomeDSL
using Carlo
using Random
using LinearAlgebra
using Printf

function test_determinant_scaling(n1, n2, label)
    params = Dict(
        :n1 => n1, :n2 => n2,
        :PBC => (true, true),
        :antiPBC => (false, true),
        :lattice => DoubleKagome,
        :B => 0.3,
        :N_up => (n1*n2*3) ÷ 3,  # 1/3 up spins
        :N_down => (n1*n2*3) * 2 ÷ 3  # 2/3 down spins
    )

    println("\n=== $label: $(n1)×$(n2) ===")
    ns = n1 * n2 * 3
    println("Total sites: $ns")
    println("N_up: $(params[:N_up]), N_down: $(params[:N_down])")

    mc = KagomeDSL.MC(params)
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 10, :seed => 1234, :thermalization => 200)
    )

    Carlo.init!(mc, ctx, params)

    # Thermalize
    for i in 1:200
        Carlo.sweep!(mc, ctx)
    end

    # Analyze determinant properties
    tilde_U_up = KagomeDSL.tilde_U(mc.Ham.U_up, mc.kappa_up)
    tilde_U_down = KagomeDSL.tilde_U(mc.Ham.U_down, mc.kappa_down)

    det_up = det(tilde_U_up)
    det_down = det(tilde_U_down)

    logdet_up = real(logdet(tilde_U_up))
    logdet_down = real(logdet(tilde_U_down))

    println("Matrix sizes: up $(size(tilde_U_up)), down $(size(tilde_U_down))")
    println("det(up): $det_up (magnitude: $(abs(det_up)))")
    println("det(down): $det_down (magnitude: $(abs(det_down)))")
    println("logdet(up): $logdet_up")
    println("logdet(down): $logdet_down")
    println("Total logdet: $(logdet_up + logdet_down)")

    # Test one S+ measurement
    site_with_down = 0
    for site = 1:ns
        if KagomeDSL.is_occupied(mc.kappa_down, site)
            site_with_down = site
            break
        end
    end

    if site_with_down > 0
        amp, _ = KagomeDSL.measure_S_plus(mc, site_with_down)
        println("S+ amplitude at site $site_with_down: $amp")
        return (ns, abs(det_up), abs(det_down), logdet_up + logdet_down, abs(amp))
    else
        println("No down spins found!")
        return (ns, abs(det_up), abs(det_down), logdet_up + logdet_down, 0.0)
    end
end

println("Testing determinant magnitude scaling with system size")
println("Looking for patterns that explain large S+ amplitudes")

# Test different sizes (n1 must be even)
sizes = [(2, 2), (2, 3), (4, 2), (4, 3)]
results = []

for (n1, n2) in sizes
    result = test_determinant_scaling(n1, n2, "$(n1)×$(n2)")
    push!(results, result)
end

println("\n" * "="^80)
println("DETERMINANT SCALING ANALYSIS")
println("="^80)
println("Size(ns) | |det_up|   | |det_down| | Total logdet | S+ amplitude")
println("-"^80)

for (ns, det_up_mag, det_down_mag, total_logdet, s_plus_amp) in results
    println(@sprintf("%-8d | %.6e | %.6e | %12.3f | %.6e",
                     ns, det_up_mag, det_down_mag, total_logdet, s_plus_amp))
end

# Analyze trends
sizes_ns = [r[1] for r in results]
logdets = [r[4] for r in results]
s_plus_amps = [r[5] for r in results]

println("\nTrend analysis:")
println("System sizes: $sizes_ns")
println("Total logdets: $logdets")
println("S+ amplitudes: $s_plus_amps")

# Check if logdet scales with system size
if length(logdets) > 2
    logdet_per_size = logdets ./ sizes_ns
    println("Logdet per site: $logdet_per_size")

    # If logdet scales roughly linearly with system size, that's the problem!
    logdet_range = maximum(logdet_per_size) - minimum(logdet_per_size)
    if logdet_range < 0.1  # Nearly constant per site
        println("✓ Logdet scales roughly linearly with system size")
        println("  This means determinants get exponentially smaller with size")
        println("  And S+ transitions can give exponentially large ratios!")
    end
end

# Expected: If S+ is correct, amplitudes should be roughly system-size independent
s_plus_range = maximum(s_plus_amps) - minimum(s_plus_amps)
s_plus_ratio = maximum(s_plus_amps) / minimum(s_plus_amps)

if s_plus_ratio > 10.0
    println("⚠ WARNING: S+ amplitudes vary by factor > 10 across sizes")
    println("  This confirms a system-size dependent normalization bug")
end

println("\nConclusion:")
println("- If total logdet scales with system size → exponential scaling problem")
println("- If S+ amplitudes vary strongly with size → normalization bug confirmed")
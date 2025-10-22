#!/usr/bin/env julia --project
# Test with large system like job.jl but optimal BC

using KagomeDSL
using Carlo
using Random
using Statistics

# Larger system with optimal boundary conditions
params = Dict(
    :n1 => 4,
    :n2 => 4,
    :PBC => (true, true),
    :antiPBC => (true, false),  # optimal BC for energy
    :N_up => 24,
    :N_down => 24
)

println("Testing S+ with large system (4×4):")
println("Total sites: $(params[:n1] * params[:n2] * 3)")
println("Half-filling: N_up=$(params[:N_up]), N_down=$(params[:N_down])")

mc = KagomeDSL.MC(params)
ctx = Carlo.MCContext{Random.Xoshiro}(
    Dict(:binsize => 20, :seed => 1234, :thermalization => 2000)
)

Carlo.init!(mc, ctx, params)

# Longer thermalization for larger system
println("Thermalizing...")
for i in 1:2000
    Carlo.sweep!(mc, ctx)
end

# Collect measurements
println("Measuring...")
s_plus_amps = Float64[]

for i in 1:1000
    Carlo.sweep!(mc, ctx)

    if i % 10 == 0
        ns = params[:n1] * params[:n2] * 3
        s_plus_amp_total = 0.0 + 0.0im

        for site = 1:ns
            amp, _ = KagomeDSL.measure_S_plus(mc, site)
            s_plus_amp_total += amp
        end

        push!(s_plus_amps, real(s_plus_amp_total / ns))
    end
end

s_plus_mean = mean(s_plus_amps)
s_plus_err = std(s_plus_amps) / sqrt(length(s_plus_amps))

println("\n=== LARGE SYSTEM RESULTS ===")
println("S+ amplitude: $s_plus_mean ± $s_plus_err")
println("Expected: ~0.667")
println("Ratio to expected: $(s_plus_mean / 0.667)")

# Try different seeds to check consistency
println("\nChecking with different random seed...")
Random.seed!(5678)
ctx2 = Carlo.MCContext{Random.Xoshiro}(
    Dict(:binsize => 20, :seed => 5678, :thermalization => 1000)
)
mc2 = KagomeDSL.MC(params)
Carlo.init!(mc2, ctx2, params)

for i in 1:1000
    Carlo.sweep!(mc2, ctx2)
end

# Quick measurement
ns = params[:n1] * params[:n2] * 3
s_plus_total = 0.0 + 0.0im
for site = 1:ns
    amp, _ = KagomeDSL.measure_S_plus(mc2, site)
    s_plus_total += amp
end
s_plus_alt = real(s_plus_total / ns)

println("Alternative seed result: $s_plus_alt")
println("Consistent: $(abs(s_plus_mean - s_plus_alt) < 0.1)")
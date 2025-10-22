#!/usr/bin/env julia --project
# Test with exact job.jl parameters

using KagomeDSL
using Carlo
using Random
using Statistics

# Exact parameters from job.jl
params = Dict(
    :n1 => 4,
    :n2 => 4,
    :PBC => (true, true),
    :antiPBC => (false, true),  # job.jl setting
    :N_up => 24,  # ns ÷ 2 where ns = 4*4*3 = 48
    :N_down => 24
)

println("Testing with exact job.jl parameters:")
println("n1=$(params[:n1]), n2=$(params[:n2])")
println("PBC=$(params[:PBC]), antiPBC=$(params[:antiPBC])")
println("N_up=$(params[:N_up]), N_down=$(params[:N_down])")
println("Total sites: $(params[:n1] * params[:n2] * 3)")

mc = KagomeDSL.MC(params)
ctx = Carlo.MCContext{Random.Xoshiro}(
    Dict(:binsize => 10, :seed => 1234, :thermalization => 1000)
)

Carlo.init!(mc, ctx, params)

# Longer thermalization like job.jl
println("Thermalizing for 1000 sweeps...")
for i in 1:1000
    Carlo.sweep!(mc, ctx)
end

# Measure energy
println("Measuring energy...")
energies = Float64[]
for i in 1:2000
    Carlo.sweep!(mc, ctx)
    if i % 20 == 0
        OL = KagomeDSL.getOL(mc, mc.kappa_up, mc.kappa_down)
        ns = params[:n1] * params[:n2] * 3
        push!(energies, OL / ns)
    end
end

energy_per_site = mean(energies)
error = std(energies) / sqrt(length(energies))

println("\n=== FINAL RESULT (job.jl parameters) ===")
println("Energy per site: $energy_per_site ± $error")
println("Expected: -0.4286")
println("Difference: $(abs(energy_per_site + 0.4286))")

# Check if S+ measurements still work in this larger system
println("\nTesting S+ measurements...")
ns = params[:n1] * params[:n2] * 3
s_plus_amps = ComplexF64[]
for i in 1:10  # Sample a few measurements
    Carlo.sweep!(mc, ctx)
    s_plus_amp = 0.0 + 0.0im
    for site = 1:ns
        amp, _ = KagomeDSL.measure_S_plus(mc, site)
        s_plus_amp += amp
    end
    push!(s_plus_amps, s_plus_amp / ns)
end

println("S+ amplitude (avg): $(mean(real.(s_plus_amps))) ± $(std(real.(s_plus_amps))/sqrt(length(s_plus_amps)))")
println("S+ measurements working: $(all(isfinite.(s_plus_amps)))")
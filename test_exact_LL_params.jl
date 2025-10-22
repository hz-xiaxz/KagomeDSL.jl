#!/usr/bin/env julia --project
# Test with exact LL parameters: B field, balanced, correct BC

using KagomeDSL
using Carlo
using Random
using Statistics

# Exact LL parameters from CLAUDE.md
n1, n2 = 4, 4  # 4n rule
ns = n1 * n2 * 3  # 48 sites  
N_up = N_down = 24  # Balanced

# LL boundary conditions
PBC = (true, true)
antiPBC = (false, true)

# Try without B field first (B=0 for ν=1/2 Laughlin state)
B = 0.0

params = Dict(
    :n1 => n1,
    :n2 => n2,
    :PBC => PBC,
    :antiPBC => antiPBC,
    :N_up => N_up,
    :N_down => N_down
)

println("=== Testing with exact LL boundary conditions ===")
println("System: $(n1)×$(n2) with N_up=$N_up, N_down=$N_down (balanced)")
println("PBC=$PBC, antiPBC=$antiPBC")
println("B field = $B")

mc = KagomeDSL.MC(params)
ctx = Carlo.MCContext{Random.Xoshiro}(
    Dict(:binsize => 10, :seed => 1234, :thermalization => 3000)
)

Carlo.init!(mc, ctx, params)

# Long thermalization
println("Thermalizing...")
for i in 1:3000
    Carlo.sweep!(mc, ctx)
end

# Measurements
println("Collecting measurements...")
energies = Float64[]
s_plus_measures = Dict(
    "norm_ndown" => Float64[],
    "avg_sq" => Float64[]
)

for i in 1:5000
    Carlo.sweep!(mc, ctx)

    if i % 40 == 0
        # Energy
        OL = KagomeDSL.getOL(mc, mc.kappa_up, mc.kappa_down)
        push!(energies, real(OL) / ns)

        # S+ measurements
        s_plus_sq_sum = 0.0

        for site = 1:ns
            amp, amp_sq = KagomeDSL.measure_S_plus(mc, site)
            s_plus_sq_sum += amp_sq
        end

        push!(s_plus_measures["avg_sq"], s_plus_sq_sum / ns)
        push!(s_plus_measures["norm_ndown"], s_plus_sq_sum / N_down)
    end
end

# Results
println("\n=== RESULTS ===")
e_mean = mean(energies)
e_err = std(energies) / sqrt(length(energies))
println("Energy per site: $e_mean ± $e_err")

println()
val_ndown = mean(s_plus_measures["norm_ndown"])
err_ndown = std(s_plus_measures["norm_ndown"]) / sqrt(length(s_plus_measures["norm_ndown"]))
println("(sum ⟨|S+|²⟩) / N_down = $val_ndown ± $err_ndown")

val_avg = mean(s_plus_measures["avg_sq"])
err_avg = std(s_plus_measures["avg_sq"]) / sqrt(length(s_plus_measures["avg_sq"]))
println("(sum ⟨|S+|²⟩) / ns = $val_avg ± $err_avg")

println()
println("Target: 0.667 ≈ 2/3")
println("Difference from 2/3: $(abs(val_ndown - 2/3))")

#!/usr/bin/env julia --project
# Test S+ in Landau Level (LL) regime with B field

using KagomeDSL
using Carlo
using Random
using Statistics

# For LL physics, need B ≠ 0 and balanced N_up = N_down
n1, n2 = 4, 4  # 4n rule
ns = n1 * n2 * 3  # 48 sites
N_up = N_down = 24  # Balanced, half filling

# Calculate B for LL: B = imbalance * π / (n1 * n2 * 2√3)
# For zero imbalance (balanced), let's try a small field
imbalance = 0
B = imbalance * π / (n1 * n2 * 2 * sqrt(3))

# Or use a reasonable field strength
B = 0.1  # Non-zero field for LL physics

params = Dict(
    :n1 => n1,
    :n2 => n2,
    :PBC => (true, true),
    :antiPBC => (false, true),  # LL boundary conditions
    :N_up => N_up,
    :N_down => N_down,
    :B => B
)

println("=== Landau Level S+ Analysis ===")
println("System: $(n1)×$(n2) with N_up=$N_up, N_down=$N_down")
println("Magnetic field B = $B")
println("PBC=$(params[:PBC]), antiPBC=$(params[:antiPBC])")

mc = KagomeDSL.MC(params)
ctx = Carlo.MCContext{Random.Xoshiro}(
    Dict(:binsize => 10, :seed => 1234, :thermalization => 2000)
)

Carlo.init!(mc, ctx, params)

# Longer thermalization for larger system
println("Thermalizing...")
for i in 1:2000
    Carlo.sweep!(mc, ctx)
end

# Collect measurements
println("Collecting measurements...")
energies = Float64[]
s_plus_measures = Dict(
    "avg" => Float64[],
    "avg_sq" => Float64[],
    "norm_ndown" => Float64[],
    "sum" => Float64[],
    "sum_sq" => Float64[]
)

for i in 1:3000
    Carlo.sweep!(mc, ctx)

    if i % 30 == 0
        # Energy
        OL = KagomeDSL.getOL(mc, mc.kappa_up, mc.kappa_down)
        push!(energies, OL / ns)

        # S+ measurements
        s_plus_sum = 0.0
        s_plus_sq_sum = 0.0

        for site = 1:ns
            amp, amp_sq = KagomeDSL.measure_S_plus(mc, site)
            s_plus_sum += abs(amp)
            s_plus_sq_sum += amp_sq
        end

        push!(s_plus_measures["avg"], s_plus_sum / ns)
        push!(s_plus_measures["avg_sq"], s_plus_sq_sum / ns)
        push!(s_plus_measures["norm_ndown"], s_plus_sq_sum / N_down)
        push!(s_plus_measures["sum"], s_plus_sum)
        push!(s_plus_measures["sum_sq"], s_plus_sq_sum)
    end
end

# Results
println("\n=== RESULTS ===")
println("Energy per site: $(mean(energies)) ± $(std(energies)/sqrt(length(energies)))")
println()
println("S+ measures:")
for (key, label) in [
    ("avg", "(sum |⟨S+⟩|) / ns"),
    ("avg_sq", "(sum ⟨|S+|²⟩) / ns"),
    ("norm_ndown", "(sum ⟨|S+|²⟩) / N_down"),
    ("sum", "sum |⟨S+⟩|"),
    ("sum_sq", "sum ⟨|S+|²⟩")
]
    val = mean(s_plus_measures[key])
    err = std(s_plus_measures[key]) / sqrt(length(s_plus_measures[key]))
    diff_from_2_3 = abs(val - 2/3)
    println("$label = $val ± $err (diff from 2/3: $diff_from_2_3)")
end

println("\nExpected: ~0.667 ≈ 2/3")

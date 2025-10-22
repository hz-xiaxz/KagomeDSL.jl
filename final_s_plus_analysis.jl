#!/usr/bin/env julia --project
# Comprehensive S+ analysis to find the 0.667 value

using KagomeDSL
using Carlo
using Random
using Statistics

# Use optimal boundary conditions
params = Dict(
    :n1 => 4,
    :n2 => 2,
    :PBC => (true, true),
    :antiPBC => (true, false),
    :N_up => 12,
    :N_down => 12
)

println("=== Comprehensive S+ Analysis ===")
println("System: $(params[:n1])×$(params[:n2]) with N_up=$(params[:N_up]), N_down=$(params[:N_down])")

mc = KagomeDSL.MC(params)
ctx = Carlo.MCContext{Random.Xoshiro}(
    Dict(:binsize => 10, :seed => 1234, :thermalization => 1000)
)

Carlo.init!(mc, ctx, params)

# Thermalize
println("Thermalizing...")
for i in 1:1000
    Carlo.sweep!(mc, ctx)
end

# Collect detailed measurements
println("Collecting measurements...")
ns = params[:n1] * params[:n2] * 3
N_down = params[:N_down]

# Track different quantities
energies = Float64[]
sum_s_plus = Float64[]        # sum_i |⟨S+⟩_i|
sum_s_plus_sq = Float64[]     # sum_i ⟨|S+|²⟩_i
avg_s_plus = Float64[]        # (sum_i |⟨S+⟩_i|) / ns
avg_s_plus_sq = Float64[]     # (sum_i ⟨|S+|²⟩_i) / ns
norm_by_ndown = Float64[]     # (sum_i ⟨|S+|²⟩_i) / N_down

for i in 1:2000
    Carlo.sweep!(mc, ctx)

    if i % 20 == 0
        # Energy
        OL = KagomeDSL.getOL(mc, mc.kappa_up, mc.kappa_down)
        push!(energies, OL / ns)

        # S+ measurements
        s_plus_amp_sum = 0.0
        s_plus_sq_sum = 0.0

        for site = 1:ns
            amp, amp_sq = KagomeDSL.measure_S_plus(mc, site)
            s_plus_amp_sum += abs(amp)
            s_plus_sq_sum += amp_sq
        end

        push!(sum_s_plus, s_plus_amp_sum)
        push!(sum_s_plus_sq, s_plus_sq_sum)
        push!(avg_s_plus, s_plus_amp_sum / ns)
        push!(avg_s_plus_sq, s_plus_sq_sum / ns)
        push!(norm_by_ndown, s_plus_sq_sum / N_down)
    end
end

# Print all results
println("\n=== RESULTS ===")
println("Energy per site: $(mean(energies)) ± $(std(energies)/sqrt(length(energies)))")
println()
println("Different S+ measures:")
println("1. (sum_i |⟨S+⟩_i|) / ns = $(mean(avg_s_plus)) ± $(std(avg_s_plus)/sqrt(length(avg_s_plus)))")
println("2. (sum_i ⟨|S+|²⟩_i) / ns = $(mean(avg_s_plus_sq)) ± $(std(avg_s_plus_sq)/sqrt(length(avg_s_plus_sq)))")
println("3. (sum_i ⟨|S+|²⟩_i) / N_down = $(mean(norm_by_ndown)) ± $(std(norm_by_ndown)/sqrt(length(norm_by_ndown)))")
println("4. sum_i |⟨S+⟩_i| = $(mean(sum_s_plus)) ± $(std(sum_s_plus)/sqrt(length(sum_s_plus)))")
println("5. sum_i ⟨|S+|²⟩_i = $(mean(sum_s_plus_sq)) ± $(std(sum_s_plus_sq)/sqrt(length(sum_s_plus_sq)))")
println()
println("Expected value: ~0.667 ≈ 2/3")
println()
println("Which one is closest to 2/3?")
for (i, (label, val)) in enumerate([
    ("(sum |⟨S+⟩|) / ns", mean(avg_s_plus)),
    ("(sum ⟨|S+|²⟩) / ns", mean(avg_s_plus_sq)),
    ("(sum ⟨|S+|²⟩) / N_down", mean(norm_by_ndown)),
    ("sum |⟨S+⟩|", mean(sum_s_plus)),
    ("sum ⟨|S+|²⟩", mean(sum_s_plus_sq))
])
    println("$i. $label = $val (diff from 2/3 = $(abs(val - 2/3)))")
end

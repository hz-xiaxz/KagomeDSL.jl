#!/usr/bin/env julia --project
# Validation of S+ measurements showing the 2/3 sum rule
# 
# This script demonstrates that with correct Landau Level (LL) parameters:
# - Balanced system (N_up = N_down)
# - LL boundary conditions: PBC=(true,true), antiPBC=(false,true)
# - Zero magnetic field (B=0) or zero imbalance
# - System size following 4n rule
#
# The S+ operator measurements satisfy the sum rule:
#   (sum_i ⟨|S+_i|²⟩) / N_down ≈ 2/3 ≈ 0.667

using KagomeDSL
using Carlo
using Random
using Statistics

println("="^60)
println("Validating S+ measurement: sum rule = 2/3")
println("="^60)

# Landau Level parameters
n1, n2 = 4, 4  # Must follow 4n rule: n1 = n2 = 4*n
ns = n1 * n2 * 3  # Total sites = 48
N_up = N_down = 24  # Balanced: N_up = N_down

params = Dict(
    :n1 => n1,
    :n2 => n2,
    :PBC => (true, true),      # Periodic in both directions
    :antiPBC => (false, true),  # Antiperiodic only in y (LL condition)
    :N_up => N_up,
    :N_down => N_down,
    :B => 0.0  # Zero magnetic field
)

println("\nSystem parameters:")
println("  Lattice size: $(n1)×$(n2) = $ns sites")
println("  Spinons: N_up=$N_up, N_down=$N_down (balanced)")
println("  Boundary: PBC=$(params[:PBC]), antiPBC=$(params[:antiPBC])")
println("  Magnetic field: B=$(params[:B])")

# Initialize Monte Carlo
mc = KagomeDSL.MC(params)
ctx = Carlo.MCContext{Random.Xoshiro}(
    Dict(:binsize => 10, :seed => 1234, :thermalization => 3000)
)
Carlo.init!(mc, ctx, params)

# Thermalization
println("\nThermalizing (3000 sweeps)...")
for i in 1:3000
    Carlo.sweep!(mc, ctx)
end

# Measurements
println("Collecting measurements (5000 sweeps)...")
energies = Float64[]
s_plus_sum_rule = Float64[]  # (sum ⟨|S+|²⟩) / N_down

n_measure = 5000
measure_interval = 40

for i in 1:n_measure
    Carlo.sweep!(mc, ctx)

    if i % measure_interval == 0
        # Energy measurement
        OL = KagomeDSL.getOL(mc, mc.kappa_up, mc.kappa_down)
        push!(energies, OL / ns)

        # S+ measurement: sum over all sites
        s_plus_sq_total = 0.0
        for site = 1:ns
            amp, amp_sq = KagomeDSL.measure_S_plus(mc, site)
            s_plus_sq_total += amp_sq
        end

        # Normalize by number of down spinons
        push!(s_plus_sum_rule, s_plus_sq_total / N_down)
    end
end

# Results
println("\n" * "="^60)
println("RESULTS")
println("="^60)

e_mean = mean(energies)
e_err = std(energies) / sqrt(length(energies))
println("\nEnergy per site:")
println("  Measured: $e_mean ± $e_err")
println("  Expected: -0.4286 (Heisenberg on Kagome)")

s_mean = mean(s_plus_sum_rule)
s_err = std(s_plus_sum_rule) / sqrt(length(s_plus_sum_rule))
diff_from_2_3 = abs(s_mean - 2/3)

println("\nS+ sum rule: (∑ᵢ ⟨|S⁺ᵢ|²⟩) / N_down")
println("  Measured: $s_mean ± $s_err")
println("  Expected: 0.6667 (2/3)")
println("  Difference: $diff_from_2_3")

if diff_from_2_3 < 0.01
    println("\n✓ SUCCESS: Sum rule validated! (within 1% of 2/3)")
else
    println("\n⚠ WARNING: Deviation from 2/3 is $(diff_from_2_3*100)%")
    println("  Consider: longer thermalization, more measurements, or larger system")
end

println("\n" * "="^60)

#!/usr/bin/env julia --project
# Test S+ measurements with correct boundary conditions

using KagomeDSL
using Carlo
using Random
using Statistics

# Use the boundary conditions that give closest energy to -0.4286
params = Dict(
    :n1 => 4,
    :n2 => 2,
    :PBC => (true, true),
    :antiPBC => (true, false),  # antiPBC in x-direction - closest to expected energy
    :N_up => 12,
    :N_down => 12
)

println("Testing S+ measurements with optimal boundary conditions:")
println("PBC=$(params[:PBC]), antiPBC=$(params[:antiPBC])")
println("System: $(params[:n1])×$(params[:n2]) with N_up=$(params[:N_up]), N_down=$(params[:N_down])")

mc = KagomeDSL.MC(params)
ctx = Carlo.MCContext{Random.Xoshiro}(
    Dict(:binsize => 10, :seed => 1234, :thermalization => 1000)
)

Carlo.init!(mc, ctx, params)

# Thermalize well
println("Thermalizing...")
for i in 1:1000
    Carlo.sweep!(mc, ctx)
end

# Collect measurements
println("Collecting measurements...")
energies = Float64[]
s_plus_amps = Float64[]
s_plus_amps_sq = Float64[]

for i in 1:2000
    Carlo.sweep!(mc, ctx)

    if i % 20 == 0
        # Energy measurement
        OL = KagomeDSL.getOL(mc, mc.kappa_up, mc.kappa_down)
        ns = params[:n1] * params[:n2] * 3
        push!(energies, OL / ns)

        # S+ measurements
        s_plus_amp_total = 0.0 + 0.0im
        s_plus_amp_sq_total = 0.0

        for site = 1:ns
            amp, amp_sq = KagomeDSL.measure_S_plus(mc, site)
            s_plus_amp_total += amp
            s_plus_amp_sq_total += amp_sq
        end

        push!(s_plus_amps, real(s_plus_amp_total / ns))
        push!(s_plus_amps_sq, s_plus_amp_sq_total / ns)
    end
end

# Results
energy_mean = mean(energies)
energy_err = std(energies) / sqrt(length(energies))

s_plus_mean = mean(s_plus_amps)
s_plus_err = std(s_plus_amps) / sqrt(length(s_plus_amps))

s_plus_sq_mean = mean(s_plus_amps_sq)
s_plus_sq_err = std(s_plus_amps_sq) / sqrt(length(s_plus_amps_sq))

println("\n=== RESULTS WITH OPTIMAL BOUNDARY CONDITIONS ===")
println("Energy per site: $energy_mean ± $energy_err")
println("Expected energy: -0.4286")
println()
println("S+ amplitude: $s_plus_mean ± $s_plus_err")
println("S+ amplitude²: $s_plus_sq_mean ± $s_plus_sq_err")
println("Expected S+ amplitude: ~0.667")
println()
println("Ratio S+²/|S+|²: $(s_plus_sq_mean / s_plus_mean^2)")

# Check final configuration
println("\nFinal configuration:")
n_up_occupied = count(!iszero, mc.kappa_up)
n_down_occupied = count(!iszero, mc.kappa_down)
println("Up spinons: $n_up_occupied sites occupied")
println("Down spinons: $n_down_occupied sites occupied")
println("Total: $(n_up_occupied + n_down_occupied) / $(length(mc.kappa_up)) sites")
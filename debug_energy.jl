#!/usr/bin/env julia --project
# Test energy calculation with different parameters

using KagomeDSL
using Carlo
using Random
using Statistics

function test_energy(params_dict, label)
    println("\n=== Testing: $label ===")
    println("Parameters: $(params_dict)")

    mc = KagomeDSL.MC(params_dict)
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 10, :seed => 1234, :thermalization => 200)
    )

    Carlo.init!(mc, ctx, params_dict)

    # Thermalize
    for i in 1:200
        Carlo.sweep!(mc, ctx)
    end

    # Measure
    energies = Float64[]
    for i in 1:500
        Carlo.sweep!(mc, ctx)
        if i % 10 == 0
            OL = KagomeDSL.getOL(mc, mc.kappa_up, mc.kappa_down)
            ns = params_dict[:n1] * params_dict[:n2] * 3
            push!(energies, OL / ns)
        end
    end

    energy_per_site = mean(energies)
    error = std(energies) / sqrt(length(energies))

    println("Energy per site: $energy_per_site ± $error")
    println("kappa_up: $(mc.kappa_up)")
    println("kappa_down: $(mc.kappa_down)")

    return energy_per_site
end

# Test 1: Small system (what we ran before)
params1 = Dict(
    :n1 => 2, :n2 => 1, :PBC => (false, false),
    :N_up => 3, :N_down => 3
)

# Test 2: Larger system with PBC (might be closer to expected)
params2 = Dict(
    :n1 => 4, :n2 => 2, :PBC => (true, true),
    :N_up => 12, :N_down => 12
)

# Test 3: Different boundary conditions
params3 = Dict(
    :n1 => 2, :n2 => 2, :PBC => (true, true),
    :N_up => 6, :N_down => 6
)

energy1 = test_energy(params1, "Small 2×1 open BC")
energy2 = test_energy(params2, "Large 4×2 PBC")
energy3 = test_energy(params3, "Medium 2×2 PBC")

println("\n=== SUMMARY ===")
println("2×1 open: $energy1")
println("4×2 PBC:  $energy2")
println("2×2 PBC:  $energy3")
println("Expected: -0.4286")
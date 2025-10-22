#!/usr/bin/env julia --project
# Test energy with PBC and antiPBC combinations

using KagomeDSL
using Carlo
using Random
using Statistics

function test_energy_bc(params_dict, label)
    println("\n=== Testing: $label ===")
    println("PBC: $(params_dict[:PBC])")
    println("antiPBC: $(params_dict[:antiPBC])")

    mc = KagomeDSL.MC(params_dict)
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 10, :seed => 1234, :thermalization => 500)
    )

    Carlo.init!(mc, ctx, params_dict)

    # Thermalize longer
    for i in 1:500
        Carlo.sweep!(mc, ctx)
    end

    # Measure
    energies = Float64[]
    for i in 1:1000
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
    return energy_per_site
end

# Test different boundary condition combinations
# System size that should give good thermodynamic limit

# Test 1: PBC in both directions
params1 = Dict(
    :n1 => 4, :n2 => 2,
    :PBC => (true, true),
    :antiPBC => (false, false),
    :N_up => 12, :N_down => 12
)

# Test 2: PBC with antiPBC in one direction
params2 = Dict(
    :n1 => 4, :n2 => 2,
    :PBC => (true, true),
    :antiPBC => (false, true),  # antiPBC in y direction
    :N_up => 12, :N_down => 12
)

# Test 3: PBC with antiPBC in other direction
params3 = Dict(
    :n1 => 4, :n2 => 2,
    :PBC => (true, true),
    :antiPBC => (true, false),  # antiPBC in x direction
    :N_up => 12, :N_down => 12
)

# Test 4: PBC with antiPBC in both directions
params4 = Dict(
    :n1 => 4, :n2 => 2,
    :PBC => (true, true),
    :antiPBC => (true, true),   # antiPBC in both directions
    :N_up => 12, :N_down => 12
)

energy1 = test_energy_bc(params1, "PBC only")
energy2 = test_energy_bc(params2, "PBC + antiPBC(y)")
energy3 = test_energy_bc(params3, "PBC + antiPBC(x)")
energy4 = test_energy_bc(params4, "PBC + antiPBC(both)")

println("\n=== BOUNDARY CONDITION COMPARISON ===")
println("PBC only:          $energy1")
println("PBC + antiPBC(y):  $energy2")
println("PBC + antiPBC(x):  $energy3")
println("PBC + antiPBC(both): $energy4")
println("Expected:          -0.4286")

# Find which one matches best
differences = [
    abs(energy1 + 0.4286),
    abs(energy2 + 0.4286),
    abs(energy3 + 0.4286),
    abs(energy4 + 0.4286)
]
best_idx = argmin(differences)
labels = ["PBC only", "PBC + antiPBC(y)", "PBC + antiPBC(x)", "PBC + antiPBC(both)"]
println("\nClosest match: $(labels[best_idx]) with difference $(differences[best_idx])")
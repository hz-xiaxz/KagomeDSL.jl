#!/usr/bin/env julia --project
# Debug S+ measurement implementation

using KagomeDSL
using Carlo
using Random
using Statistics

# Test with simple system first
params = Dict(
    :n1 => 2, :n2 => 1,
    :PBC => (false, false),
    :N_up => 2, :N_down => 4  # Asymmetric to see different S+ behavior
)

println("Debugging S+ measurement:")
println("System: $(params[:n1])×$(params[:n2]) with N_up=$(params[:N_up]), N_down=$(params[:N_down])")

mc = KagomeDSL.MC(params)
ctx = Carlo.MCContext{Random.Xoshiro}(
    Dict(:binsize => 10, :seed => 1234, :thermalization => 100)
)

Carlo.init!(mc, ctx, params)

# Short thermalization
for i in 1:100
    Carlo.sweep!(mc, ctx)
end

println("\nConfiguration after thermalization:")
println("kappa_up:   $(mc.kappa_up)")
println("kappa_down: $(mc.kappa_down)")

# Test S+ measurement on individual sites
ns = length(mc.kappa_up)
println("\nS+ measurements by site:")
for site = 1:ns
    has_down = KagomeDSL.is_occupied(mc.kappa_down, site)
    amp, amp_sq = KagomeDSL.measure_S_plus(mc, site)

    println("Site $site: has_down=$has_down, S+ amp=$amp, |S+|²=$amp_sq")
end

# Average S+ amplitude
s_plus_total = 0.0 + 0.0im
s_plus_sq_total = 0.0

for site = 1:ns
    amp, amp_sq = KagomeDSL.measure_S_plus(mc, site)
    global s_plus_total += amp
    global s_plus_sq_total += amp_sq
end

avg_s_plus = s_plus_total / ns
avg_s_plus_sq = s_plus_sq_total / ns

println("\nAveraged over all sites:")
println("⟨S+⟩ = $(real(avg_s_plus))")
println("⟨|S+|²⟩ = $avg_s_plus_sq")

# Test different N_up/N_down ratios
println("\n=== Testing different particle number ratios ===")

test_cases = [
    (2, 4),  # More down spins
    (3, 3),  # Equal
    (4, 2),  # More up spins
]

for (N_up, N_down) in test_cases
    params_test = Dict(
        :n1 => 2, :n2 => 1,
        :PBC => (false, false),
        :N_up => N_up, :N_down => N_down
    )

    mc_test = KagomeDSL.MC(params_test)
    ctx_test = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 10, :seed => 1234, :thermalization => 50)
    )

    Carlo.init!(mc_test, ctx_test, params_test)

    for i in 1:50
        Carlo.sweep!(mc_test, ctx_test)
    end

    # Measure S+
    ns_test = length(mc_test.kappa_up)
    s_plus_sum = 0.0 + 0.0im

    for site = 1:ns_test
        amp, _ = KagomeDSL.measure_S_plus(mc_test, site)
        s_plus_sum += amp
    end

    avg_amp = real(s_plus_sum / ns_test)
    println("N_up=$N_up, N_down=$N_down: ⟨S+⟩ = $avg_amp")
end
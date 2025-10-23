#!/usr/bin/env julia --project
# Benchmark S+ measurement performance

using KagomeDSL
using Carlo
using Random
using BenchmarkTools
using LinearAlgebra

println("="^60)
println("Benchmarking S+ Measurement Performance")
println("="^60)

# Setup MC state
params = Dict(
    :n1 => 4, :n2 => 4,
    :PBC => (true, true),
    :antiPBC => (false, true),
    :N_up => 24, :N_down => 24
)

mc = KagomeDSL.MC(params)
ctx = Carlo.MCContext{Random.Xoshiro}(
    Dict(:binsize => 10, :seed => 1234, :thermalization => 100)
)
Carlo.init!(mc, ctx, params)

# Thermalize
for i in 1:100
    Carlo.sweep!(mc, ctx)
end

println("\nSystem: $(params[:n1])×$(params[:n2]), N_up=$(params[:N_up]), N_down=$(params[:N_down])")
println("Total sites: $(params[:n1] * params[:n2] * 3)")

# Benchmark individual S+ measurement
println("\n--- Single Site S+ Measurement ---")
site = 10
@btime KagomeDSL.measure_S_plus($mc, $site)

# Benchmark full lattice S+ measurement
println("\n--- Full Lattice S+ Sum ---")
function measure_full_s_plus(mc)
    ns = length(mc.kappa_up)
    s_plus_sq_sum = 0.0
    for site = 1:ns
        _, amp_sq = KagomeDSL.measure_S_plus(mc, site)
        s_plus_sq_sum += amp_sq
    end
    return s_plus_sq_sum
end

@btime measure_full_s_plus($mc)

# Benchmark tilde_U construction (current bottleneck)
println("\n--- tilde_U Matrix Construction ---")
@btime KagomeDSL.tilde_U($mc.Ham.U_up, $mc.kappa_up)

# Benchmark logdet calculation
println("\n--- logdet Calculation ---")
tilde_U_test = KagomeDSL.tilde_U(mc.Ham.U_up, mc.kappa_up)
@btime real(logdet($tilde_U_test))

println("\n" * "="^60)
println("Performance Summary:")
println("- Main bottleneck: tilde_U construction in get_log_det_ratio")
println("- Called 4 times per S+ measurement (2x for mc_n, 2x for mc_np1)")
println("- Optimization target: Cache tilde_U or log_det values")
println("="^60)

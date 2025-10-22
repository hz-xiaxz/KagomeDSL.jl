#!/usr/bin/env julia --project
# Quick test script to see S+ measurement results

using KagomeDSL
using Carlo
using Random
using Statistics

# Small system for quick testing
params = Dict(
    :n1 => 2,
    :n2 => 1,
    :PBC => (false, false),
    :N_up => 3,
    :N_down => 3,
    :seed => 1234,
    :thermalization => 100,
    :sweeps => 500,
    :binsize => 10
)

println("Setting up Monte Carlo simulation...")
println("System: $(params[:n1])×$(params[:n2]) lattice with N_up=$(params[:N_up]), N_down=$(params[:N_down])")

# Create MC state
mc = KagomeDSL.MC(params)

# Set up Carlo context
ctx = Carlo.MCContext{Random.Xoshiro}(
    Dict(
        :binsize => params[:binsize],
        :seed => params[:seed],
        :thermalization => params[:thermalization]
    )
)

# Initialize
println("Initializing...")
Carlo.init!(mc, ctx, params)

# Run thermalization
println("Thermalizing...")
for i in 1:params[:thermalization]
    Carlo.sweep!(mc, ctx)
end

# Run measurements
println("Running measurements...")
measurements = Dict()
measurements[:S_plus_amp] = Float64[]
measurements[:S_plus_amp_sq] = Float64[]
measurements[:energy] = Float64[]

for i in 1:params[:sweeps]
    Carlo.sweep!(mc, ctx)

    # Measure every 10 sweeps
    if i % 10 == 0
        # Measure energy
        OL = KagomeDSL.getOL(mc, mc.kappa_up, mc.kappa_down)
        push!(measurements[:energy], OL / (params[:n1] * params[:n2] * 3))

        # Measure S+ amplitude
        ns = length(mc.kappa_up)
        s_plus_amp = 0.0 + 0.0im
        s_plus_amp_sq = 0.0
        for site = 1:ns
            amp, amp_sq = KagomeDSL.measure_S_plus(mc, site)
            s_plus_amp += amp
            s_plus_amp_sq += amp_sq
        end
        push!(measurements[:S_plus_amp], real(s_plus_amp / ns))
        push!(measurements[:S_plus_amp_sq], s_plus_amp_sq / ns)
    end
end

# Analyze results
println("\n=== RESULTS ===")
println("Energy per site: $(mean(measurements[:energy])) ± $(std(measurements[:energy]) / sqrt(length(measurements[:energy])))")
println("S+ amplitude: $(mean(measurements[:S_plus_amp])) ± $(std(measurements[:S_plus_amp]) / sqrt(length(measurements[:S_plus_amp])))")
println("S+ amplitude squared: $(mean(measurements[:S_plus_amp_sq])) ± $(std(measurements[:S_plus_amp_sq]) / sqrt(length(measurements[:S_plus_amp_sq])))")

println("\nConfiguration:")
println("kappa_up: $(mc.kappa_up)")
println("kappa_down: $(mc.kappa_down)")
#!/usr/bin/env julia --project
# Analyze S+ measurements more carefully

using KagomeDSL
using Carlo
using Random
using Statistics

# Test with optimal boundary conditions
params = Dict(
    :n1 => 4, :n2 => 2,
    :PBC => (true, true),
    :antiPBC => (true, false),  # optimal boundary conditions
    :N_up => 12, :N_down => 12
)

println("Analyzing S+ measurements with optimal boundary conditions:")
println("Half-filling case: N_up=$(params[:N_up]), N_down=$(params[:N_down])")

mc = KagomeDSL.MC(params)
ctx = Carlo.MCContext{Random.Xoshiro}(
    Dict(:binsize => 10, :seed => 1234, :thermalization => 1000)
)

Carlo.init!(mc, ctx, params)

# Thermalize well
for i in 1:1000
    Carlo.sweep!(mc, ctx)
end

# Analyze S+ measurements
println("\nAnalyzing configuration:")
println("kappa_up:   $(mc.kappa_up)")
println("kappa_down: $(mc.kappa_down)")

ns = params[:n1] * params[:n2] * 3
n_down_occupied = count(!iszero, mc.kappa_down)
n_up_occupied = count(!iszero, mc.kappa_up)

println("Sites with down spins: $n_down_occupied / $ns")
println("Sites with up spins: $n_up_occupied / $ns")

# Measure S+ at all sites
all_amps = ComplexF64[]
sites_with_down = Int[]
amps_with_down = ComplexF64[]

for site = 1:ns
    has_down = KagomeDSL.is_occupied(mc.kappa_down, site)
    amp, _ = KagomeDSL.measure_S_plus(mc, site)

    push!(all_amps, amp)
    if has_down
        push!(sites_with_down, site)
        push!(amps_with_down, amp)
    end
end

# Statistics
avg_all_sites = real(mean(all_amps))
avg_sites_with_down = length(amps_with_down) > 0 ? real(mean(amps_with_down)) : 0.0

println("\n=== S+ MEASUREMENT ANALYSIS ===")
println("Average over ALL sites: $avg_all_sites")
println("Average over sites WITH down spins: $avg_sites_with_down")
println("Number of sites with down spins: $(length(sites_with_down))")
println("Expected individual amplitude: ~0.667")

# Show individual measurements for sites with down spins
if length(sites_with_down) > 0
    println("\nIndividual S+ amplitudes at sites with down spins:")
    for (i, site) in enumerate(sites_with_down[1:min(10, end)])
        amp = real(amps_with_down[i])
        println("  Site $site: $amp")
    end
end

# Test different configurations by sweeping
println("\n=== TESTING MULTIPLE CONFIGURATIONS ===")
all_site_averages = Float64[]
down_site_averages = Float64[]

for sweep_i in 1:20
    Carlo.sweep!(mc, ctx)

    # Measure S+ for this configuration
    all_amps_sweep = ComplexF64[]
    down_amps_sweep = ComplexF64[]

    for site = 1:ns
        has_down = KagomeDSL.is_occupied(mc.kappa_down, site)
        amp, _ = KagomeDSL.measure_S_plus(mc, site)

        push!(all_amps_sweep, amp)
        if has_down
            push!(down_amps_sweep, amp)
        end
    end

    push!(all_site_averages, real(mean(all_amps_sweep)))
    if length(down_amps_sweep) > 0
        push!(down_site_averages, real(mean(down_amps_sweep)))
    end
end

println("Average S+ (all sites) over 20 sweeps: $(mean(all_site_averages)) ± $(std(all_site_averages))")
if length(down_site_averages) > 0
    println("Average S+ (down sites) over 20 sweeps: $(mean(down_site_averages)) ± $(std(down_site_averages))")
end
#!/usr/bin/env julia --project
# Debug the determinant ratio calculation to find normalization issue

using KagomeDSL
using Carlo
using Random
using LinearAlgebra

# Test with small system
params = Dict(
    :n1 => 2, :n2 => 2,
    :PBC => (true, true),
    :antiPBC => (false, true),
    :lattice => DoubleKagome,
    :B => 0.5,
    :N_up => 4, :N_down => 8  # More down spins to ensure we find them
)

println("Debugging S+ amplitude calculation")
println("System: $(params[:n1])×$(params[:n2]), ns = $(params[:n1]*params[:n2]*3)")
println("N_up = $(params[:N_up]), N_down = $(params[:N_down])")

mc = KagomeDSL.MC(params)
ctx = Carlo.MCContext{Random.Xoshiro}(
    Dict(:binsize => 10, :seed => 1234, :thermalization => 100)
)

Carlo.init!(mc, ctx, params)

# Thermalize briefly
for i in 1:100
    Carlo.sweep!(mc, ctx)
end

# Find a site with down spin
global site_with_down = 0
for site = 1:length(mc.kappa_down)
    if KagomeDSL.is_occupied(mc.kappa_down, site)
        global site_with_down = site
        break
    end
end

if site_with_down == 0
    println("No down spins found!")
    exit(1)
end

println("\nAnalyzing S+ transition at site $site_with_down")
println("Initial configuration:")
println("  kappa_up:   $(mc.kappa_up)")
println("  kappa_down: $(mc.kappa_down)")

# Perform the transition manually to debug
mc_np1 = KagomeDSL.spin_plus_transition(mc, site_with_down)

println("\nAfter S+ transition:")
println("  kappa_up:   $(mc_np1.kappa_up)")
println("  kappa_down: $(mc_np1.kappa_down)")

# Check the determinant components
println("\nDebugging determinant calculation:")

# Original state determinants
tilde_U_up_n = KagomeDSL.tilde_U(mc.Ham.U_up, mc.kappa_up)
tilde_U_down_n = KagomeDSL.tilde_U(mc.Ham.U_down, mc.kappa_down)

println("Original state (n):")
println("  tilde_U_up size: $(size(tilde_U_up_n))")
println("  tilde_U_down size: $(size(tilde_U_down_n))")

det_up_n = det(tilde_U_up_n)
det_down_n = det(tilde_U_down_n)
logdet_up_n = real(logdet(tilde_U_up_n))
logdet_down_n = real(logdet(tilde_U_down_n))

println("  det(tilde_U_up): $det_up_n")
println("  det(tilde_U_down): $det_down_n")
println("  logdet(tilde_U_up): $logdet_up_n")
println("  logdet(tilde_U_down): $logdet_down_n")

# Transition state determinants
tilde_U_up_np1 = KagomeDSL.tilde_U(mc_np1.Ham.U_up, mc_np1.kappa_up)
tilde_U_down_np1 = KagomeDSL.tilde_U(mc_np1.Ham.U_down, mc_np1.kappa_down)

println("\nTransition state (n+1):")
println("  tilde_U_up size: $(size(tilde_U_up_np1))")
println("  tilde_U_down size: $(size(tilde_U_down_np1))")

det_up_np1 = det(tilde_U_up_np1)
det_down_np1 = det(tilde_U_down_np1)
logdet_up_np1 = real(logdet(tilde_U_up_np1))
logdet_down_np1 = real(logdet(tilde_U_down_np1))

println("  det(tilde_U_up): $det_up_np1")
println("  det(tilde_U_down): $det_down_np1")
println("  logdet(tilde_U_up): $logdet_up_np1")
println("  logdet(tilde_U_down): $logdet_down_np1")

# Calculate the ratio
log_ratio = KagomeDSL.get_log_det_ratio(mc, mc_np1)
ratio = exp(log_ratio)

println("\nRatio calculation:")
println("  log_det_n = $(logdet_up_n + logdet_down_n)")
println("  log_det_np1 = $(logdet_up_np1 + logdet_down_np1)")
println("  log_ratio = $log_ratio")
println("  ratio = exp(log_ratio) = $ratio")

# Compare with direct determinant ratio
det_ratio_direct = (det_up_np1 * det_down_np1) / (det_up_n * det_down_n)
println("  Direct det ratio = $det_ratio_direct")

# Current S+ amplitude
current_amplitude, current_amp_sq = KagomeDSL.measure_S_plus(mc, site_with_down)
println("\nCurrent S+ amplitude: $current_amplitude")
println("Expected range: ~0.1 to ~2.0")

if abs(current_amplitude) > 5.0
    println("⚠ WARNING: Amplitude suspiciously large!")
end

# Test different normalization approaches
println("\nTesting potential normalizations:")

# 1. Divide by system size
normalized_1 = current_amplitude / (params[:n1] * params[:n2] * 3)
println("  1. Amplitude / ns = $normalized_1")

# 2. Take square root (since det ratio might be squared)
normalized_2 = sqrt(abs(current_amplitude)) * sign(real(current_amplitude))
println("  2. sqrt(|amplitude|) = $normalized_2")

# 3. Normalize by number of particles
N_total = params[:N_up] + params[:N_down]
normalized_3 = current_amplitude / N_total
println("  3. Amplitude / N_total = $normalized_3")

# 4. Logarithmic scaling correction
normalized_4 = current_amplitude / exp(N_total / 10.0)
println("  4. Amplitude / exp(N/10) = $normalized_4")
using KagomeDSL
using Random
using Statistics
using LinearAlgebra

Random.seed!(1234)

# Test with small system - 2×1 Kagome lattice
n1, n2 = 2, 1
ns = n1 * n2 * 3  # 6 sites
N_up, N_down = 2, 4
J = 1.0
B = 0.0
PBC = (true, true)
antiPBC = (true, false)

lattice = kagome_lattice(n1, n2; PBC=PBC, antiPBC=antiPBC)
Ham = Hamiltonian{N_up,N_down}(lattice, J, B)
mc = MCState(Ham, lattice)

# Thermalize
n_therm = 1000
for _ in 1:n_therm
    sweep!(mc)
end

# Measure S+ in different ways
n_measure = 5000
s_plus_site = zeros(ComplexF64, ns)
s_plus_sq_site = zeros(Float64, ns)

for _ in 1:n_measure
    sweep!(mc)

    for site in 1:ns
        s_plus_amp, s_plus_sq = measure_S_plus(mc, site)
        s_plus_site[site] += s_plus_amp
        s_plus_sq_site[site] += s_plus_sq
    end
end

s_plus_site ./= n_measure
s_plus_sq_site ./= n_measure

println("\n=== Analysis of S+ measurements ===")
println("System: $(n1)×$(n2) with N_up=$N_up, N_down=$N_down")
println("\nBy site averages:")
for site in 1:ns
    println("Site $site: ⟨S+⟩=$(abs(s_plus_site[site])), ⟨|S+|²⟩=$(s_plus_sq_site[site])")
end

# Different averaging schemes
avg_all_sites = mean(abs.(s_plus_site))
sum_all_sites = sum(abs.(s_plus_site))
sum_sq_all_sites = sum(s_plus_sq_site)

println("\n=== Different averaging schemes ===")
println("1. Average over all sites: $(avg_all_sites)")
println("2. Sum over all sites: $(sum_all_sites)")
println("3. Sum over all sites / N_down: $(sum_all_sites / N_down)")
println("4. Sum of ⟨|S+|²⟩ over all sites: $(sum_sq_all_sites)")
println("5. Sum of ⟨|S+|²⟩ / N_down: $(sum_sq_all_sites / N_down)")

# Check if any combination gives ~0.667
println("\n=== Checking for 0.667 ≈ 2/3 ===")
println("Is sum_sq / N_down ≈ 2/3? $(sum_sq_all_sites / N_down) (expected 0.667)")

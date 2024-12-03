using Random

@testset "Hamiltonian" begin
    DK = DoubleKagome(1.0, 4, 3, (false, false))
    H = KagomeDSL.Hmat(DK)
    # check H is Hermitian
    @test H ≈ H'
    @test isempty(findall(x -> !(x ≈ 0), H - H'))
    # test that the first row of H has only 2,3 position elements =1 and others =0
    @test H[1, 1] == 0
    @test H[1, 2] == -1
    @test H[1, 3] == -1
    for i = 4:36
        @test H[1, i] == 0
    end
    @test H[3, 13] == 1

    DK2 = KagomeDSL.DoubleKagome(1.0, 4, 3, (true, false))
    H2 = KagomeDSL.Hmat(DK2)
    @test isempty(findall(x -> !(x ≈ 0), H2 - H2'))
    @test H2 ≈ H2'

    @test H2[1, 2] == -1
    @test H2[1, 3] == -1
    @test CartesianIndex(1, 11) in DK2.nn
    @test H2[1, 11] == -1 # horizontal PBC
    @test H2[1, 3*4+1] == 0 # vertical OBC

    # other boundaries should be tested

end

@testset "orbitals" begin
    DK2 = DoubleKagome(1.0, 4, 3, (true, false))
    χ = 1.0
    N_up = 18
    N_down = 18
    Han = KagomeDSL.Hamiltonian(N_up, N_down, DK2)
    U_up = Han.U_up
    U_down = Han.U_down
    num = ns(DK2)
    @test size(U_up) == (num, N_up)
    @test size(U_down) == (num, N_down)
end

@testset "getxprime" begin
    # TODO More careful tests here
    DK = DoubleKagome(1.0, 4, 3, (false, false))
    N_up = 1
    N_down = 0
    ham = KagomeDSL.Hamiltonian(N_up, N_down, DK)
    kappa_up = vcat([1], fill(0, 35))
    kappa_down = vcat([0], collect(1:35))
    # Sx Sy will flip 1 to spin down, the nearest neighbor will be spin up
    # this gives 2 configurations
    xprime = KagomeDSL.getxprime(ham, kappa_up, kappa_down)
    @test length(keys(xprime)) == 3
    # Sz interaction
    Sz_sum = (length(DK.nn) - 2) * (1 / 2)
    Sz_sum += 2 * (2 * (-1 / 4))
    @test xprime[(-1, -1, -1, -1)] == Sz_sum
    # Sx Sy interaction, only happens in the bonds having site 1
    Sx_Sy_sum = 0.0
    # Sx Sy = 1/2()
    # x1 is flip 1 down, filp 2 up
    kappa_up1 = vcat([0], [1], fill(0, 34))
    kappa_down1 = vcat([1], [0], collect(2:35))
    k_up = 2
    l_up = 1
    k_down = 1
    l_down = 1
    conf = (k_up, l_up, k_down, l_down)
    @test conf in keys(xprime)
    @test xprime[conf] == 1 / 2 * 2
    # x2 is flip 1 up, flip 3 down
    k_up_3 = 3
    l_up_3 = 1
    k_down_3 = 1
    l_down_3 = 2
    conf_3 = (k_up_3, l_up_3, k_down_3, l_down_3)
    @test conf_3 in keys(xprime)
    @test xprime[conf_3] == 1 / 2 * 2
end

@testset "getOL" begin
    # consider up: [1,0,1,1,1,1,1,0,0,0,0,0] down: [1,0,1,1,1,1,1,0,0,0,0,0]
    kappa_up = vcat([1], [0], collect(2:6), fill(0, 5))
    mc = KagomeDSL.MC(
        Dict(:n1 => 2, :n2 => 2, :PBC => (false, false), :N_up => 6, :N_down => 6),
    )
    # consider up: [1,0,1,1,1,1,1,0,0,0,0,0] down: [0,1,0,0,0,0,0,1,1,1,1,1]
    kappa_down = vcat([0], [1], fill(0, 5), collect(2:6))
    @test KagomeDSL.getOL(mc, kappa_up, kappa_down) != 0.0
end

@testset "Sz tests" begin
    @testset "Basic spin states" begin
        kappa_up = [1, 0, 2]
        kappa_down = [0, 2, 0]

        # Test up spin
        @test Sz(1, kappa_up, kappa_down) == 0.5

        # Test down spin
        @test Sz(2, kappa_up, kappa_down) == -0.5
    end

    @testset "Error conditions" begin
        kappa_up = [1, 2, 0]
        kappa_down = [0, 2, 1]

        # Test double occupation
        @test_throws ArgumentError Sz(2, kappa_up, kappa_down)

        # Test bounds error
        @test_throws BoundsError Sz(4, kappa_up, kappa_down)
        @test_throws BoundsError Sz(0, kappa_up, kappa_down)

        # Test mismatched vector lengths
        kappa_up_long = [1, 2, 0, 1]
        @test_throws DimensionMismatch Sz(1, kappa_up_long, kappa_down)
    end

    @testset "Edge cases" begin
        # Single site case
        @test Sz(1, [1], [0]) == 0.5
        @test Sz(1, [0], [1]) == -0.5
        @test_throws ArgumentError Sz(1, [0], [0])
        @test_throws ArgumentError Sz(1, [1], [1])
    end

    @testset "Type stability" begin
        kappa_up = [1, 0]
        kappa_down = [0, 1]

        # Test that return type is always Float64
        @inferred Sz(1, kappa_up, kappa_down)
        @inferred Sz(2, kappa_up, kappa_down)
    end

    @testset "Large system" begin
        n = 1000
        kappa_up = zeros(Int, n)
        kappa_down = zeros(Int, n)

        # Set some spins
        kappa_up[1] = 1     # up spin at start
        kappa_down[end] = 1 # down spin at end

        @test Sz(1, kappa_up, kappa_down) == 0.5
        @test Sz(n, kappa_up, kappa_down) == -0.5
        @test_throws ArgumentError Sz(2, kappa_up, kappa_down) # empty site
    end

    @testset "Random configurations" begin
        n = 100
        rng = Random.MersenneTwister(42)

        # Create random valid configuration
        sites = shuffle(rng, 1:n)
        half_n = n ÷ 2

        kappa_up = zeros(Int, n)
        kappa_down = zeros(Int, n)

        # Assign first half to up spins, second half to down spins
        for (i, site) in enumerate(sites[1:half_n])
            kappa_up[site] = i
        end
        for (i, site) in enumerate(sites[half_n+1:end])
            kappa_down[site] = i
        end

        # Test random sites
        for site in sites
            if !iszero(kappa_up[site])
                @test Sz(site, kappa_up, kappa_down) == 0.5
            else
                @test Sz(site, kappa_up, kappa_down) == -0.5
            end
        end
    end
end

@testset "spinInteraction! tests" begin
    @testset "Basic spin flips" begin
        # Initialize test configuration
        kappa_up = [1, 0, 2]   # up spins at sites 1 and 3
        kappa_down = [0, 1, 0]  # down spin at site 2
        xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()

        # Test S+_i S-_j: site 2(↓) -> site 1(↑)
        spinInteraction!(xprime, kappa_up, kappa_down, 2, 1)
        @test haskey(xprime, (2, 1, 1, 1))  # new configuration
        @test xprime[(2, 1, 1, 1)] ≈ 0.5    # coefficient should be 1/2

        # Test S-_i S+_j: site 1(↑) -> site 2(↓)
        xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()
        spinInteraction!(xprime, kappa_up, kappa_down, 1, 2)
        @test haskey(xprime, (2, 1, 1, 1))  # new configuration
        @test xprime[(2, 1, 1, 1)] ≈ 0.5    # coefficient should be 1/2
    end

    @testset "No action cases" begin
        kappa_up = [1, 0, 2]
        kappa_down = [0, 1, 0]
        xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()

        # Test when both sites have same spin
        spinInteraction!(xprime, kappa_up, kappa_down, 1, 3)  # both up
        @test isempty(xprime)

        spinInteraction!(xprime, kappa_up, kappa_down, 2, 3)
        @test !isempty(xprime)
        @test xprime[(2, 2, 3, 1)] == 0.5
    end

    @testset "Multiple interactions" begin
        kappa_up = [1, 0, 2]
        kappa_down = [0, 1, 0]
        xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()

        # Apply same interaction twice
        spinInteraction!(xprime, kappa_up, kappa_down, 1, 2)
        spinInteraction!(xprime, kappa_up, kappa_down, 1, 2)
        @test xprime[(2, 1, 1, 1)] ≈ 1.0  # coefficients should add
    end

    @testset "Edge cases" begin
        # Single site case
        kappa_up = [1]
        kappa_down = [0]
        xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()
        spinInteraction!(xprime, kappa_up, kappa_down, 1, 1)
        @test isempty(xprime)  # no self-interaction

        # Empty configuration
        kappa_up = Int[]
        kappa_down = Int[]
        xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()
        @test_throws BoundsError spinInteraction!(xprime, kappa_up, kappa_down, 1, 1)
    end

    @testset "Large system" begin
        n = 1000
        kappa_up = zeros(Int, n)
        kappa_down = zeros(Int, n)
        kappa_up[1] = 1    # up spin at first site
        kappa_down[2] = 1  # down spin at second site

        xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()
        spinInteraction!(xprime, kappa_up, kappa_down, 1, 2)

        @test haskey(xprime, (2, 1, 1, 1))
        @test xprime[(2, 1, 1, 1)] ≈ 0.5
    end

    @testset "Accumulation behavior" begin
        kappa_up = [1, 0, 2]
        kappa_down = [0, 1, 0]
        xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()

        # Add some initial value
        xprime[(2, 1, 1, 1)] = 0.25

        # Apply interaction
        spinInteraction!(xprime, kappa_up, kappa_down, 1, 2)
        @test xprime[(2, 1, 1, 1)] ≈ 0.75  # 0.25 + 0.5
    end

    @testset "Type stability" begin
        kappa_up = [1, 0]
        kappa_down = [0, 1]
        xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()

        # Test that the function maintains type stability
        @inferred spinInteraction!(xprime, kappa_up, kappa_down, 1, 2)
    end
end

@testset "DoubleKagome Boundary Conditions" begin
    @testset "get_boundary_shifts" begin
        n1, n2 = 4, 3  # n1 must be even for DoubleKagome
        ns = n1 * n2 * 3
        lat = DoubleKagome(1.0, n1, n2, (true, true); antiPBC = (true, false))
        tunneling = zeros(Float64, ns, ns)

        # Example link dictionaries
        link_in = Dict((1, 2) => 1.0, (2, 1) => 1.0)

        link_inter = Dict((1, -1) => 1.0, (-1, 1) => 1.0)

        # Test in-cell tunneling (should not be affected by antiperiodic BC)
        s1, s2 = 1, 2
        apply_boundary_conditions!(tunneling, lat, s1, s2, n1, ns, link_in, link_inter)
        @test tunneling[s1, s2] ≈ 1.0
        @test tunneling[s2, s1] ≈ 1.0

        # Test inter-cell tunneling with antiperiodic BC
        tunneling = zeros(Float64, ns, ns)
        s1, s2 = 1, 4
        apply_boundary_conditions!(tunneling, lat, s1, s2, n1, ns, link_in, link_inter)

        # Check if crossing boundary in first direction gives negative sign
        boundary_cross = abs(s1 - s2) > 3n1
        expected_sign = boundary_cross ? -1.0 : 1.0
        if haskey(link_inter, bondNum(s1, s2))
            @test tunneling[s1, s2] ≈ expected_sign * link_inter[bondNum(s1, s2)]
        end
    end
end

@testset "bondNum Function" begin
    # Test cases for positive indices
    @testset "Positive Indices" begin
        @test bondNum(1, 7) == (1, 7)
        @test bondNum(6, 12) == (6, 12)
        @test bondNum(7, 13) == (1, 7)
        @test bondNum(12, 18) == (6, 12)
        @test bondNum(13, 19) == (1, 7)
    end

    # Test cases for negative indices
    @testset "Negative Indices" begin
        @test bondNum(-1, 1) == (5, 7)
        @test bondNum(1, -1) == (1, -1)
        @test bondNum(-6, 0) == (6, 12)
        @test bondNum(-7, -1) == (5, 11)
        @test bondNum(-12, -6) == (6, 12)
        @test bondNum(-13, -7) == (5, 11)
    end

    # Test cases for edge cases
    @testset "Edge Cases" begin
        @test bondNum(0, 6) == (6, 12)
        @test bondNum(-1, 0) == (5, 6)
        @test bondNum(1, 0) == (1, 0)
        @test bondNum(6, 0) == (6, 0)
    end
end

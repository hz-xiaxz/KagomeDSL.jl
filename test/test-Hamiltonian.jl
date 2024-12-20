using Random
using KagomeDSL
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

    # Anti-periodic boundary condition tests
    @testset "Anti-periodic boundary conditions" begin
        # Test anti-PBC in x direction
        DK_antiX = DoubleKagome(1.0, 4, 3, (true, false); antiPBC = (true, false))
        H_antiX = KagomeDSL.Hmat(DK_antiX)
        @test H_antiX ≈ H_antiX'  # Check Hermiticity
        @test H_antiX[1, 2] == -1  # Normal in-cell bond
        @test H_antiX[1, 3] == -1  # Normal in-cell bond
        @test H_antiX[5, 7] == -1   # Normal inter-cell bond
        @test H_antiX[1, 11] == 1 # Boundary crossing with anti-PBC in x

        # Test anti-PBC in y direction
        DK_antiY = DoubleKagome(1.0, 4, 3, (true, true); antiPBC = (false, true))
        H_antiY = KagomeDSL.Hmat(DK_antiY)
        @test H_antiY ≈ H_antiY'  # Check Hermiticity
        @test H_antiY[1, 11] == -1   # Normal horizontal bond
        @test H_antiY[1, 3+4*6] == -1  # Vertical bond with anti-PBC

        # Test anti-PBC in both directions
        DK_antiBoth = DoubleKagome(1.0, 4, 3, (true, true); antiPBC = (true, true))
        H_antiBoth = KagomeDSL.Hmat(DK_antiBoth)
        @test H_antiBoth ≈ H_antiBoth'  # Check Hermiticity
        @test H_antiBoth[1, 11] == 1  # x-direction anti-PBC
        @test H_antiBoth[1, 27] == -1  # y-direction anti-PBC
        @test H_antiBoth[11, 27] == 1  # Double crossing should give positive sign

        # Test mixed boundary conditions
        DK_mixed = DoubleKagome(1.0, 4, 3, (true, true); antiPBC = (true, false))
        H_mixed = KagomeDSL.Hmat(DK_mixed)
        @test H_mixed ≈ H_mixed'  # Check Hermiticity
        @test H_mixed[1, 11] == 1  # x-direction anti-PBC
        @test H_mixed[1, 3+4*6] == 1  # y-direction normal PBC
    end

    # Test error cases
    @testset "Error cases" begin
        # Test invalid anti-PBC configuration
        @test_throws ArgumentError DoubleKagome(
            1.0,
            4,
            3,
            (false, false);
            antiPBC = (true, false),
        )
        @test_throws ArgumentError DoubleKagome(
            1.0,
            4,
            3,
            (false, false);
            antiPBC = (false, true),
        )
    end
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
    @test xprime[conf] == -1 / 2 * 2
    # x2 is flip 1 up, flip 3 down
    k_up_3 = 3
    l_up_3 = 1
    k_down_3 = 1
    l_down_3 = 2
    conf_3 = (k_up_3, l_up_3, k_down_3, l_down_3)
    @test conf_3 in keys(xprime)
    @test xprime[conf_3] == -1 / 2 * 2
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
        @test xprime[(2, 1, 1, 1)] ≈ -0.5    # coefficient should be -1/2

        # Test S-_i S+_j: site 1(↑) -> site 2(↓)
        xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()
        spinInteraction!(xprime, kappa_up, kappa_down, 1, 2)
        @test haskey(xprime, (2, 1, 1, 1))  # new configuration
        @test xprime[(2, 1, 1, 1)] ≈ -0.5    # coefficient should be -1/2
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
        @test xprime[(2, 2, 3, 1)] == -0.5
    end

    @testset "Multiple interactions" begin
        kappa_up = [1, 0, 2]
        kappa_down = [0, 1, 0]
        xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()

        # Apply same interaction twice
        spinInteraction!(xprime, kappa_up, kappa_down, 1, 2)
        spinInteraction!(xprime, kappa_up, kappa_down, 1, 2)
        @test xprime[(2, 1, 1, 1)] ≈ -1.0  # coefficients should add
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
        @test xprime[(2, 1, 1, 1)] ≈ -0.5
    end

    @testset "Accumulation behavior" begin
        kappa_up = [1, 0, 2]
        kappa_down = [0, 1, 0]
        xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()

        # Add some initial value
        xprime[(2, 1, 1, 1)] = 0.25

        # Apply interaction
        spinInteraction!(xprime, kappa_up, kappa_down, 1, 2)
        @test xprime[(2, 1, 1, 1)] ≈ -0.25  # 0.25 - 0.5
    end

    @testset "Type stability" begin
        kappa_up = [1, 0]
        kappa_down = [0, 1]
        xprime = Dict{Tuple{Int,Int,Int,Int},Float64}()

        # Test that the function maintains type stability
        @inferred spinInteraction!(xprime, kappa_up, kappa_down, 1, 2)
    end
end
@testset "apply_boundary_conditions!" begin
    @testset "Basic boundary conditions" begin
        n1, n2 = 4, 3  # n1 must be even for DoubleKagome
        ns = n1 * n2 * 6  # 6 sites per unit cell

        # Define test cases with different boundary conditions
        test_cases = [
            (
                "PBC in both directions",
                (true, true),
                (false, false),
                Dict((1, 5, -1, 0) => 1.0, (5, 1, 1, 0) => 1.0),
                [(1, 11, 1.0), (11, 1, 1.0)],  # (s1, s2, expected_value)
            ),
            (
                "Anti-PBC in x direction",
                (true, true),
                (true, false),
                Dict((1, 5, -1, 0) => 1.0, (5, 1, 1, 0) => 1.0),
                [(1, 11, -1.0), (11, 1, -1.0)],
            ),
            (
                "Anti-PBC in y direction",
                (true, true),
                (false, true),
                Dict((3, 6, 0, -1) => 1.0, (6, 3, 0, 1) => 1.0),
                [(3, 3 * 2 * n1 + 6, -1.0), (3 * 2 * n1 + 6, 3, -1.0)],
            ),
        ]

        for (name, PBC, antiPBC, link_inter, tests) in test_cases
            @testset "$name" begin
                lat = DoubleKagome(1.0, n1, n2, PBC; antiPBC = antiPBC)

                for (s1, s2, expected) in tests
                    tunneling = zeros(Float64, ns, ns)
                    apply_boundary_conditions!(tunneling, lat, s1, s2, link_inter)
                    @test tunneling[s1, s2] ≈ expected
                end
            end
        end
    end

    @testset "Multiple bonds and crossings" begin
        n1, n2 = 4, 3
        ns = n1 * n2 * 3
        lat = DoubleKagome(1.0, n1, n2, (true, true); antiPBC = (true, true))

        # Test diagonal crossings
        link_inter = Dict((1, 6, -1, -1) => 1.0)
        tunneling = zeros(Float64, ns, ns)
        s1, s2 = 1, ns
        apply_boundary_conditions!(tunneling, lat, s1, s2, link_inter)
        @test tunneling[s1, s2] ≈ 1.0  # Signs cancel for diagonal crossing
    end

    @testset "Edge cases" begin
        n1, n2 = 4, 3
        ns = n1 * n2 * 3
        lat = DoubleKagome(1.0, n1, n2, (true, true); antiPBC = (true, true))


        # Test non-existent bonds
        link_inter = Dict((1, 5, 1, 0) => 1.0)
        tunneling = zeros(Float64, ns, ns)
        s1, s2 = 1, 7  # Sites that don't have a defined bond
        apply_boundary_conditions!(tunneling, lat, s1, s2, link_inter)
        @test tunneling[s1, s2] ≈ 0.0  # No bond should be added
    end

    @testset "Error handling" begin
        n1, n2 = 4, 3
        ns = n1 * n2 * 3
        lat = DoubleKagome(1.0, n1, n2, (true, true))
        tunneling = zeros(Float64, ns, ns)
        link_inter = Dict((1, 5, -1, 0) => 1.0)

        # Test same cell error
        s1, s2 = 1, 2  # Same unit cell
        @test_throws AssertionError apply_boundary_conditions!(
            tunneling,
            lat,
            s1,
            s2,
            link_inter,
        )

        # Test invalid site indices
        @test_throws AssertionError apply_boundary_conditions!(
            tunneling,
            lat,
            0,
            1,
            link_inter,
        )
        @test_throws AssertionError apply_boundary_conditions!(
            tunneling,
            lat,
            1,
            ns + 1,
            link_inter,
        )
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

        link_inter = Dict((1, 5, -1, 0) => 1.0, (5, 1, 1, 0) => 1.0)


        # Test inter-cell tunneling with antiperiodic BC
        tunneling = zeros(Float64, ns, ns)
        s1, s2 = 1, 11
        apply_boundary_conditions!(tunneling, lat, s1, s2, link_inter)

        # Check if crossing boundary in first direction gives negative sign
        boundary_cross = true
        expected_sign = boundary_cross ? -1.0 : 1.0
        if haskey(link_inter, (1, 5, -1, 0))
            @test tunneling[s1, s2] ≈ expected_sign * link_inter[(1, 5, -1, 0)]
        end
    end
end

@testset "unitcell_coord" begin
    n1, n2 = 3, 2  # Small lattice for testing

    # Test first unit cell (sites 1-6)
    @test KagomeDSL.unitcell_coord(1, n1, n2) ≈ [0.0, 0.0]
    @test KagomeDSL.unitcell_coord(6, n1, n2) ≈ [0.0, 0.0]

    # Test middle unit cell
    middle_cell = 2  # second cell in first row
    @test KagomeDSL.unitcell_coord(6 * middle_cell - 5, n1, n2) ≈ [4.0, 0.0]

    # Test last unit cell
    last_site = n1 * n2 * 6
    expected_last = [(n1 - 1) * 4.0 + (n2 - 1) * 1.0, (n2 - 1) * sqrt(3.0)]
    @test KagomeDSL.unitcell_coord(last_site, n1, n2) ≈ expected_last
    # Test row transitions
    # Last cell in first row
    @test KagomeDSL.unitcell_coord(6 * n1, n1, n2) ≈ [(n1 - 1) * 4.0, 0.0]
    # First cell in second row
    @test KagomeDSL.unitcell_coord(6 * n1 + 1, n1, n2) ≈ [1.0, sqrt(3.0)]

    # Test error cases
    @test_throws AssertionError KagomeDSL.unitcell_coord(0, n1, n2)
    @test_throws AssertionError KagomeDSL.unitcell_coord(n1 * n2 * 6 + 1, n1, n2)

    # Test coordinate system
    # Check x-direction spacing
    x1 = KagomeDSL.unitcell_coord(1, n1, n2)[1]
    x2 = KagomeDSL.unitcell_coord(7, n1, n2)[1]  # Next cell in x-direction
    @test x2 - x1 ≈ 4.0

    # Check y-direction spacing
    y1 = KagomeDSL.unitcell_coord(1, n1, n2)[2]
    y2 = KagomeDSL.unitcell_coord(6 * n1 + 1, n1, n2)[2]  # Next cell in y-direction
    @test y2 - y1 ≈ sqrt(3.0)
end

@testset "unitcell_diff with integer rounding" begin
    # Define basis vectors
    a1 = [4.0, 0.0]
    a2 = [1.0, sqrt(3.0)]

    # Test exact unit translations
    @testset "Exact unit translations" begin
        # One unit in a1 direction
        dx, dy = KagomeDSL.unitcell_diff(a1, [0.0, 0.0])
        @test dx == 1
        @test dy == 0

        # One unit in a2 direction
        dx, dy = KagomeDSL.unitcell_diff(a2, [0.0, 0.0])
        @test dx == 0
        @test dy == 1

        # One unit in each direction
        dx, dy = KagomeDSL.unitcell_diff(a1 + a2, [0.0, 0.0])
        @test dx == 1
        @test dy == 1
    end

    # Test rounding behavior
    @testset "Rounding behavior" begin
        # Slightly perturbed coordinates should round to same integers
        for ε in [-0.1, 0.1]
            # Near one unit in a1
            dx, dy = KagomeDSL.unitcell_diff(a1 + [ε, ε], [0.0, 0.0])
            @test dx == 1
            @test dy == 0

            # Near one unit in a2
            dx, dy = KagomeDSL.unitcell_diff(a2 + [ε, ε], [0.0, 0.0])
            @test dx == 0
            @test dy == 1
        end
    end

    # Test negative translations
    @testset "Negative translations" begin
        dx, dy = KagomeDSL.unitcell_diff([0.0, 0.0], a1)
        @test dx == -1
        @test dy == 0

        dx, dy = KagomeDSL.unitcell_diff([0.0, 0.0], a2)
        @test dx == 0
        @test dy == -1
    end

    # Test multiple unit cells
    @testset "Multiple unit cells" begin
        # Two units in a1
        dx, dy = KagomeDSL.unitcell_diff(2 * a1, [0.0, 0.0])
        @test dx == 2
        @test dy == 0

        # Two units in a2
        dx, dy = KagomeDSL.unitcell_diff(2 * a2, [0.0, 0.0])
        @test dx == 0
        @test dy == 2

        # Diagonal: two units in each direction
        dx, dy = KagomeDSL.unitcell_diff(2 * a1 + 2 * a2, [0.0, 0.0])
        @test dx == 2
        @test dy == 2
    end

    # Test zero difference
    @testset "Zero difference" begin
        dx, dy = KagomeDSL.unitcell_diff([0.0, 0.0], [0.0, 0.0])
        @test dx == 0
        @test dy == 0

        # Small perturbations should still give zero
        dx, dy = KagomeDSL.unitcell_diff([0.1, 0.1], [0.0, 0.0])
        @test dx == 0
        @test dy == 0
    end

    # Test with actual lattice coordinates
    @testset "Lattice coordinates" begin
        # Define some typical lattice positions
        positions = [
            ([0.0, 0.0], [4.0, 0.0], -1, 0),    # One unit left
            ([4.0, 0.0], [0.0, 0.0], 1, 0),     # One unit right
            ([0.0, 0.0], [1.0, sqrt(3.0)], 0, -1), # One unit down
            ([1.0, sqrt(3.0)], [0.0, 0.0], 0, 1),  # One unit up
            ([5.0, sqrt(3.0)], [0.0, 0.0], 1, 1),  # Diagonal
        ]

        for (coord1, coord2, expected_dx, expected_dy) in positions
            dx, dy = KagomeDSL.unitcell_diff(coord1, coord2)
            @test dx == expected_dx
            @test dy == expected_dy
        end
    end
end

@testset "get_boundary_shifts" begin
    @testset "Basic functionality" begin
        n1, n2 = 4, 3  # n1 must be even for DoubleKagome
        lat = DoubleKagome(1.0, n1, n2, (true, true))

        # Test same position (should never happen in practice due to assertions)
        s1, s2 = 3, 3
        @test_throws AssertionError get_boundary_shifts(lat, s1, s2)

        # Test nearest neighbor in same row
        s1, s2 = 3, 9  # Sites in adjacent cells
        shifts = get_boundary_shifts(lat, s1, s2)
        @test (1, 0, 1.0) in shifts

        # Test periodic wrapping in x direction
        s1, s2 = 3, 6 * (n1 ÷ 2 - 1) + 3  # First and last cell in row
        shifts = get_boundary_shifts(lat, s1, s2)
        @test (-(n1 ÷ 2 - 1), 0, 1.0) in shifts
        @test (n1 ÷ 2 + 1, 0, 1.0) in shifts
    end

    @testset "Boundary conditions" begin
        n1, n2 = 4, 3

        # Test with no PBC
        lat_no_pbc = DoubleKagome(1.0, n1, n2, (false, false))
        s1, s2 = 3, 9
        shifts = get_boundary_shifts(lat_no_pbc, s1, s2)
        @test length(shifts) == 1  # Only direct connection
        @test shifts == [(1, 0, 1.0)]

        # Test with PBC in x only
        lat_pbc_x = DoubleKagome(1.0, n1, n2, (true, false))
        shifts = get_boundary_shifts(lat_pbc_x, s1, s2)
        @test length(shifts) > 1  # Should include periodic images in x
        @test all(s[2] == 0 for s in shifts)  # No y-shifts

        # Test with anti-PBC in x
        lat_anti_x = DoubleKagome(1.0, n1, n2, (true, false); antiPBC = (true, false))
        s1, s2 = 3, 6 * (n1 ÷ 2 - 1) + 3  # First and last cell in row
        shifts = get_boundary_shifts(lat_anti_x, s1, s2)
        @test any(s[3] == -1.0 for s in shifts)  # Should have negative signs
    end

    @testset "Edge cases and errors" begin
        n1, n2 = 4, 3
        lat = DoubleKagome(1.0, n1, n2, (true, true))

        # Test out of range indices
        @test_throws AssertionError get_boundary_shifts(lat, 0, 1)
        @test_throws AssertionError get_boundary_shifts(lat, 1, n1 * n2 * 6 + 1)

        # Test diagonal shifts
        s1, s2 = 3, 6 * n1 ÷ 2 + 9  # Diagonal neighbor
        shifts = get_boundary_shifts(lat, s1, s2)
        @test any(s[1] != 0 && s[2] != 0 for s in shifts)  # Should have diagonal shifts
    end

    @testset "Sign consistency" begin
        n1, n2 = 4, 3
        lat = DoubleKagome(1.0, n1, n2, (true, true); antiPBC = (true, true))

        # Test sign consistency for equivalent paths
        s1, s2 = 3, 6 * n1 ÷ 2 + 9
        shifts = get_boundary_shifts(lat, s1, s2)

        # Group shifts by displacement
        by_displacement = Dict{Tuple{Int,Int},Vector{Float64}}()
        for (dx, dy, sign) in shifts
            if !haskey(by_displacement, (dx, dy))
                by_displacement[(dx, dy)] = Float64[]
            end
            push!(by_displacement[(dx, dy)], sign)
        end

        # Check that all signs are consistent for each displacement
        for (_, signs) in by_displacement
            @test length(unique(signs)) == 1  # All signs should be the same
        end
    end

    @testset "Doubled unit cell handling" begin
        n1, n2 = 4, 3
        lat = DoubleKagome(1.0, n1, n2, (true, true))

        # Test that shifts account for doubled unit cell
        s1, s2 = 3, 9  # Adjacent cells
        shifts = get_boundary_shifts(lat, s1, s2)
        @test (1, 0, 1.0) in shifts  # Direct connection

        # Test wrapping with doubled unit cell
        s1, s2 = 3, 6 * (n1 ÷ 2 - 1) + 3  # First and last cell
        shifts = get_boundary_shifts(lat, s1, s2)
        @test any(abs(s[1]) == n1 ÷ 2 + 1 for s in shifts)  # Should use halved n1
    end
end

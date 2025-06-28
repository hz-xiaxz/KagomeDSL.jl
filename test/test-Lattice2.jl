
using KagomeDSL
using Test

@testset "DoubleKagome2" begin
    # Test assertion for even n1
    @test_throws AssertionError DoubleKagome2(1.0, 3, 3, (false, false))

    # Test basic functionality
    DK = DoubleKagome2(1.0, 4, 3, (false, false))
    @test ns(DK) == 36
    @test DK.a2 == [-1.0, 0.5*sqrt(3.0)*2]
end

@testset "unitcell_coord2" begin
    n1, n2 = 3, 2  # Small lattice for testing
    DK = DoubleKagome2(1.0, 2*n1, n2, (false, false))
    # Test first unit cell (sites 1-6)
    @test KagomeDSL.unitcell_coord(DK, 1) ≈ [0.0, 0.0]
    @test KagomeDSL.unitcell_coord(DK, 6) ≈ [0.0, 0.0]

    # Test middle unit cell
    middle_cell = 2  # second cell in first row
    @test KagomeDSL.unitcell_coord(DK, 6 * middle_cell - 5) ≈ [4.0, 0.0]

    # Test last unit cell
    last_site = n1 * n2 * 6
    expected_last = (n1 - 1) * DK.a1 + (n2 - 1) * DK.a2
    @test KagomeDSL.unitcell_coord(DK, last_site) ≈ expected_last
    # Test row transitions
    # Last cell in first row
    @test KagomeDSL.unitcell_coord(DK, 6 * n1) ≈ (n1 - 1) * DK.a1
    # First cell in second row
    @test KagomeDSL.unitcell_coord(DK, 6 * n1 + 1) ≈ DK.a2

    # Test error cases
    @test_throws AssertionError KagomeDSL.unitcell_coord(DK, 0)
    @test_throws AssertionError KagomeDSL.unitcell_coord(DK, n1 * n2 * 6 + 1)

    # Test coordinate system
    # Check x-direction spacing
    x1 = KagomeDSL.unitcell_coord(DK, 1)[1]
    x2 = KagomeDSL.unitcell_coord(DK, 7)[1]  # Next cell in x-direction
    @test x2 - x1 ≈ 4.0

    # Check y-direction spacing
    y1 = KagomeDSL.unitcell_coord(DK, 1)[2]
    y2 = KagomeDSL.unitcell_coord(DK, 6 * n1 + 1)[2]  # Next cell in y-direction
    @test y2 - y1 ≈ DK.a2[2]
end

@testset "unitcell_diff2 with integer rounding" begin
    DK = DoubleKagome2(1.0, 4, 3, (false, false))
    a1 = DK.a1
    a2 = DK.a2

    # Test exact unit translations
    @testset "Exact unit translations" begin
        # One unit in a1 direction
        dx, dy = KagomeDSL.unitcell_diff(DK, a1, [0.0, 0.0])
        @test dx == 1
        @test dy == 0

        # One unit in a2 direction
        dx, dy = KagomeDSL.unitcell_diff(DK, a2, [0.0, 0.0])
        @test dx == 0
        @test dy == 1

        # One unit in each direction
        dx, dy = KagomeDSL.unitcell_diff(DK, a1 + a2, [0.0, 0.0])
        @test dx == 1
        @test dy == 1
    end

    # Test rounding behavior
    @testset "Rounding behavior" begin
        # Slightly perturbed coordinates should round to same integers
        for ε in [-0.1, 0.1]
            # Near one unit in a1
            dx, dy = KagomeDSL.unitcell_diff(DK, a1 + [ε, ε], [0.0, 0.0])
            @test dx == 1
            @test dy == 0

            # Near one unit in a2
            dx, dy = KagomeDSL.unitcell_diff(DK, a2 + [ε, ε], [0.0, 0.0])
            @test dx == 0
            @test dy == 1
        end
    end
end

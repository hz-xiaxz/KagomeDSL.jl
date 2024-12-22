@testset "KagomeLattice" begin
    # Test basic functionality without antiPBC
    K1 = KagomeDSL.Kagome(1.0, 3, 3, (false, false))
    @test ns(K1) == 27
    @test size(K1.distance_matrix) == (27, 27)

    # Test diagonal elements
    for i = 1:27
        @test K1.distance_matrix[i, i] == 0.0
    end

    # Test basic distances
    @test K1.distance_matrix[1, 2] ≈ 1.0
    @test K1.distance_matrix[1, 3] ≈ 1.0
    @test K1.distance_matrix[1, 4] ≈ 2.0
    @test K1.distance_matrix[1, 5] ≈ 3.0
    @test K1.distance_matrix[1, 6] ≈ √(2.5^2 + (√3 / 2)^2)
    @test K1.distance_matrix[1, 7] ≈ 4.0
    @test K1.distance_matrix[1, 10] ≈ 2.0
    @test K1.distance_matrix[1, 13] ≈ 2√3

    # Test boundary condition validation
    @test_throws ArgumentError KagomeDSL.Kagome(
        1.0,
        3,
        3,
        (false, false);
        antiPBC = (true, false),
    )
    @test_throws ArgumentError KagomeDSL.Kagome(
        1.0,
        3,
        3,
        (false, false);
        antiPBC = (false, true),
    )
    @test_throws ArgumentError KagomeDSL.Kagome(
        1.0,
        3,
        3,
        (true, false);
        antiPBC = (false, true),
    )

    # Test valid antiPBC configurations
    K2 = KagomeDSL.Kagome(1.0, 3, 3, (true, true); antiPBC = (true, false))
    @test ns(K2) == 27
    @test size(K2.distance_matrix) == (27, 27)

    # Test diagonal elements with antiPBC
    @inbounds for i = 1:27
        @test K2.distance_matrix[i, i] == 0.0
    end

    # Test basic distances with antiPBC
    @test K2.distance_matrix[1, 2] ≈ 1.0
    @test K2.distance_matrix[1, 3] ≈ 1.0
    @test K2.distance_matrix[1, 4] ≈ 2.0
    @test K2.distance_matrix[1, 5] ≈ 3.0
    @test K2.distance_matrix[1, 6] ≈ √(2.5^2 + (√3 / 2)^2)
    @test K2.distance_matrix[1, 7] ≈ 2.0
    @test K2.distance_matrix[1, 8] ≈ 1.0
    @test K2.distance_matrix[1, 9] ≈ √3

end

@testset "DoubleKagome" begin
    # Test assertion for even n1
    @test_throws AssertionError DoubleKagome(1.0, 3, 3, (false, false))

    # Test basic functionality
    DK = DoubleKagome(1.0, 4, 3, (false, false))
    @test ns(DK) == 36
    @test size(DK.distance_matrix) == (36, 36)

    # Test diagonal elements
    @inbounds for i = 1:36
        @test DK.distance_matrix[i, i] == 0.0
    end

    # Test basic distances
    @test DK.distance_matrix[1, 2] ≈ 1.0
    @test DK.distance_matrix[1, 3] ≈ 1.0
    @test DK.distance_matrix[1, 4] ≈ 2.0
    @test DK.distance_matrix[1, 5] ≈ 3.0
    @test DK.distance_matrix[1, 6] ≈ √(2.5^2 + (√3 / 2)^2)
    @test DK.distance_matrix[1, 7] ≈ 4.0
    @test DK.distance_matrix[1, 10] ≈ 6.0
    @test DK.distance_matrix[1, 13] ≈ 2.0

    DK3 = KagomeDSL.DoubleKagome(1.0, 4, 4, (true, true); antiPBC = (true, false))
    @test DK3.distance_matrix[1, 15] == 0.0
end

@testset "nearestNeighbor" begin
    K1 = KagomeDSL.Kagome(1.0, 3, 3, (true, false))
    @test ns(K1) == 27
    nn = K1.nn

    # Test CartesianIndex ordering
    for ind in nn
        @test ind[1] < ind[2]
    end

    # Test in-cell bonds
    @test CartesianIndex(1, 2) in nn
    @test CartesianIndex(1, 3) in nn
    @test CartesianIndex(2, 3) in nn
    @test !(CartesianIndex(1, 4) in nn)

    # Test with antiPBC
    K2 = KagomeDSL.Kagome(1.0, 3, 3, (true, true); antiPBC = (true, true))
    @test ns(K2) == 27
    nn2 = K2.nn

    # Test different nearest neighbor structure
    @test length(nn) != length(nn2)
    @test all(ind in nn2 for ind in nn)
end

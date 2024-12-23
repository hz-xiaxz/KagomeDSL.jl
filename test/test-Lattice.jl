@testset "DoubleKagome" begin
    # Test assertion for even n1
    @test_throws AssertionError DoubleKagome(1.0, 3, 3, (false, false))

    # Test basic functionality
    DK = DoubleKagome(1.0, 4, 3, (false, false))
    @test ns(DK) == 36
end

@testset "nearestNeighbor" begin
    # Create lattice and Hamiltonian
    K1 = KagomeDSL.DoubleKagome(1.0, 4, 4, (true, false))
    H1 = KagomeDSL.Hmat(K1)
    @test size(H1, 1) == 48# Test system size
    nn = KagomeDSL.get_nn(H1)  # Get nearest neighbors from Hamiltonian

    # Test ordering of pairs
    for (i, j) in nn
        @test i < j  # First index should be smaller
    end

    # Test in-cell bonds (convert CartesianIndex to linear indices)
    @test (1, 2) in nn
    @test (1, 3) in nn
    @test (2, 3) in nn
    @test !((1, 4) in nn)

    # Test with antiPBC
    K2 = KagomeDSL.DoubleKagome(1.0, 4, 4, (true, true); antiPBC = (true, true))
    H2 = KagomeDSL.Hmat(K2)
    nn2 = KagomeDSL.get_nn(H2)

    # Test different nearest neighbor structure
    @test length(nn) != length(nn2)  # Different number of bonds with antiPBC
    @test all(pair in nn2 for pair in nn)  # All original bonds should still exist

    # Test non-zero elements in Hamiltonian match nearest neighbors
    for (i, j) in nn
        @test H1[i, j] != 0  # Check bond exists in Hamiltonian
    end

    for (i, j) in nn2
        @test H2[i, j] != 0  # Check bond exists in antiPBC Hamiltonian
    end
end

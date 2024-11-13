using KagomeDSL
using Random
using Test
using LinearAlgebra

@testset "MC" begin
    param = Dict(:n1 => 4, :n2 => 3, :PBC => (true, false), :N_up => 18, :N_down => 18)
    test_mc = MC(param)
    @test length(filter(!iszero, test_mc.kappa_up)) == 18
    @test length(filter(!iszero, test_mc.kappa_down)) == 18
end

@testset "KagomeDSL.init_conf tests" begin
    @testset "Basic functionality" begin
        rng = Random.MersenneTwister(42)  # Fixed seed for reproducibility
        ns = 10
        N_up = 4

        kappa_up, kappa_down = KagomeDSL.init_conf(rng, ns, N_up)

        # Test output types and sizes
        @test kappa_up isa Vector{Int}
        @test kappa_down isa Vector{Int}
        @test length(kappa_up) == ns
        @test length(kappa_down) == ns
    end

    @testset "Conservation properties" begin
        rng = Random.MersenneTwister(42)
        ns = 12
        N_up = 5

        kappa_up, kappa_down = KagomeDSL.init_conf(rng, ns, N_up)

        # Test number of non-zero elements
        @test count(!iszero, kappa_up) == N_up
        @test count(!iszero, kappa_down) == ns - N_up

        # Test that values in kappa_up are 1:N_up
        @test sort(filter(!iszero, kappa_up)) == collect(1:N_up)

        # Test that values in kappa_down are 1:(ns-N_up)
        @test sort(filter(!iszero, kappa_down)) == collect(1:(ns-N_up))
    end

    @testset "Mutual exclusivity" begin
        rng = Random.MersenneTwister(42)
        ns = 8
        N_up = 3

        kappa_up, kappa_down = KagomeDSL.init_conf(rng, ns, N_up)

        # Test that each position is either in up or down, but not both
        for i = 1:ns
            @test xor(kappa_up[i] != 0, kappa_down[i] != 0)
        end
    end

    @testset "Edge cases" begin
        rng = Random.MersenneTwister(42)

        # Test with N_up = 0
        kappa_up, kappa_down = KagomeDSL.init_conf(rng, 5, 0)
        @test all(iszero, kappa_up)
        @test all(!iszero, kappa_down)
        @test sort(kappa_down) == collect(1:5)

        # Test with N_up = ns
        kappa_up, kappa_down = KagomeDSL.init_conf(rng, 5, 5)
        @test sort(kappa_up) == collect(1:5)
        @test all(iszero, kappa_down)

        # Test with ns = 1
        kappa_up, kappa_down = KagomeDSL.init_conf(rng, 1, 1)
        @test kappa_up == [1]
        @test kappa_down == [0]
    end

    @testset "Random distribution" begin
        rng = Random.MersenneTwister(42)
        ns = 100
        N_up = 40

        # Run multiple times to check randomness
        results = Vector{Tuple{Vector{Int},Vector{Int}}}()
        for _ = 1:10
            push!(results, KagomeDSL.init_conf(rng, ns, N_up))
        end

        # Test that we get different configurations
        @test length(unique(results)) > 1

        # Test that all configurations are valid
        for (kappa_up, kappa_down) in results
            @test count(!iszero, kappa_up) == N_up
            @test count(!iszero, kappa_down) == ns - N_up
            @test sort(filter(!iszero, kappa_up)) == collect(1:N_up)
            @test sort(filter(!iszero, kappa_down)) == collect(1:(ns-N_up))
        end
    end
end

@testset "tiled_U tests" begin
    @testset "Basic functionality" begin
        # Create a simple test matrix
        U = [
            1.0 2.0 3.0
            4.0 5.0 6.0
            7.0 8.0 9.0
        ]
        kappa = [2, 3, 1]

        result = tiled_U(U, kappa)

        # Test output dimensions
        @test size(result) == (3, 3)

        # Test correct row placement
        @test result[1, :] == U[3, :]  # kappa[3] = 1
        @test result[2, :] == U[1, :]  # kappa[1] = 2
        @test result[3, :] == U[2, :]  # kappa[2] = 3
    end

    @testset "Zero kappa entries" begin
        U = [
            1.0 2.0
            3.0 4.0
        ]
        kappa = [0, 0]

        @test_throws ArgumentError tiled_U(U, kappa)
    end

    @testset "Different matrix shapes" begin
        # Test with rectangular matrix
        U = [
            1.0 2.0
            3.0 4.0
            5.0 6.0
        ]
        kappa = [1, 0, 2]

        result = tiled_U(U, kappa)

        @test size(result) == (2, 2)
        @test result[1, :] == U[1, :]
        @test result[2, :] == U[3, :]
    end

    @testset "Complex numbers" begin
        U = [
            1.0+im 2.0+2im
            3.0+3im 4.0+4im
        ]
        kappa = [2, 1]

        result = tiled_U(U, kappa)

        @test result[1, :] == U[2, :]
        @test result[2, :] == U[1, :]
    end

    @testset "Edge cases" begin
        # Empty matrix
        U = Matrix{Float64}(undef, 0, 0)
        kappa = Int[]
        result = tiled_U(U, kappa)
        @test size(result) == (0, 0)

        # Single element
        U = reshape([1.0], 1, 1)
        kappa = [1]
        result = tiled_U(U, kappa)
        @test result == U
    end

    @testset "Input preservation" begin
        U = [
            1.0 2.0
            3.0 4.0
        ]
        U_original = copy(U)
        kappa = [1, 2]
        kappa_original = copy(kappa)

        result = tiled_U(U, kappa)

        # Test that inputs weren't modified
        @test U == U_original
        @test kappa == kappa_original
    end

    @testset "Invalid inputs" begin
        U = [
            1.0 2.0
            3.0 4.0
        ]

        # Test kappa with invalid indices
        @test_throws BoundsError tiled_U(U, [3, 1])  # Index 3 is out of bounds

        # Test mismatched dimensions
        @test_throws DimensionMismatch tiled_U(U, [1, 2, 3])
    end

    @testset "Type stability" begin
        U = [
            1.0 2.0
            3.0 4.0
        ]
        kappa = [1, 2]

        # Test that output type matches input type
        result = tiled_U(U, kappa)
        @test eltype(result) == eltype(U)

        # Test with different types
        U_int = [1 2; 3 4]
        result_int = tiled_U(U_int, kappa)
        @test eltype(result_int) == eltype(U_int)
    end
end

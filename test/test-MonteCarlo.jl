using KagomeDSL: update_W!, update_W_matrices!, is_occupied, update_configurations!, tilde_U, relabel_configuration!, spin_plus_transition
using Random
using Test
using LinearAlgebra
using Carlo
using HDF5

@testset "MC" begin
    param = Dict(:n1 => 4, :n2 => 3, :PBC => (true, false), :N_up => 18, :N_down => 18)
    test_mc = KagomeDSL.MC(param)
    # kappa_up and kappa_down are initialized in Carlo.init!, so we don't check their lengths here.
end

@testset "relabel_configuration!" begin
    kappa = [0, 5, 0, 2, 3]
    relabel_configuration!(kappa)
    @test kappa == [0, 1, 0, 2, 3]

    kappa = [10, 20, 30]
    relabel_configuration!(kappa)
    @test kappa == [1, 2, 3]

    kappa = [0, 0, 0]
    relabel_configuration!(kappa)
    @test kappa == [0, 0, 0]

    kappa = [1, 0, 0]
    relabel_configuration!(kappa)
    @test kappa == [1, 0, 0]
end

@testset "spin_plus_transition" begin
    # Setup initial MC state
    params = Dict(:n1 => 2, :n2 => 1, :PBC => (false, false), :N_up => 3, :N_down => 3)
    mc = MC(params)
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 3, :seed => 123, :thermalization => 10),
    )
    Carlo.init!(mc, ctx, params)
    
    # Initialize configuration
    mc.kappa_up = [1, 2, 3, 0, 0, 0]
    mc.kappa_down = [0, 0, 0, 1, 2, 3]

    # Site to flip
    site_to_flip = 4 # Has a down spin

    # Perform transition
    new_mc = spin_plus_transition(mc, site_to_flip)

    # Check new particle numbers
    @test size(new_mc.Ham.U_up, 2) == params[:N_up] + 1
    @test size(new_mc.Ham.U_down, 2) == params[:N_down] - 1

    # Check new kappa vectors
    @test new_mc.kappa_up[site_to_flip] != 0
    @test new_mc.kappa_down[site_to_flip] == 0
    @test count(!iszero, new_mc.kappa_up) == params[:N_up] + 1
    @test count(!iszero, new_mc.kappa_down) == params[:N_down] - 1

    # Check that W matrices are re-evaluated and not all zero
    @test !iszero(new_mc.W_up)
    @test !iszero(new_mc.W_down)
end

@testset "init_conf_qr! tests" begin
    @testset "Basic functionality and non-singularity" begin
        n1 = 2
        n2 = 2
        ns = n1 * n2 * 3 # 12 sites
        N_up = 6
        N_down = 6

        # Create a mock Hamiltonian for testing
        # Ensure U_up and U_down are full rank for a non-singular tilde_U to be possible
        U_up_mock = rand(ns, N_up)
        U_down_mock = rand(ns, N_down)
        mock_ham = Hamiltonian{N_up,N_down}(
            U_up_mock,
            U_down_mock,
            zeros(ns, ns),
            [],
            zeros(ns, N_up + 1),
            zeros(ns, N_down - 1),
        )

        mc =
            MC(mock_ham, zeros(Int, ns), zeros(Int, ns), zeros(ns, N_up), zeros(ns, N_down))

        # Call the QR-based initialization
        KagomeDSL.init_conf_qr!(mc, ns, N_up, N_down)

        # Verify kappa_up and kappa_down properties
        @test count(!iszero, mc.kappa_up) == N_up
        @test count(!iszero, mc.kappa_down) == N_down
        @test all(x -> x in 1:N_up, filter(!iszero, mc.kappa_up))
        @test all(x -> x in 1:N_down, filter(!iszero, mc.kappa_down))
        @test count(iszero, mc.kappa_up) + count(iszero, mc.kappa_down) == ns * 2 - N_up - N_down

        # Verify non-singularity of tilde_U matrices
        tilde_U_up = tilde_U(mc.Ham.U_up, mc.kappa_up)
        @test abs(det(tilde_U_up)) > eps(Float64)

        if N_down > 0
            tilde_U_down = tilde_U(mc.Ham.U_down, mc.kappa_down)
            @test abs(det(tilde_U_down)) > eps(Float64)
        end
    end
end

@testset "tilde_U tests" begin
    @testset "Basic functionality" begin
        # Create a simple test matrix
        U = [
            1.0 2.0 3.0
            4.0 5.0 6.0
            7.0 8.0 9.0
        ]
        kappa = [2, 3, 1]

        result = tilde_U(U, kappa)

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

        @test_throws ArgumentError tilde_U(U, kappa)
    end

    @testset "Different matrix shapes" begin
        # Test with rectangular matrix
        U = [
            1.0 2.0
            3.0 4.0
            5.0 6.0
        ]
        kappa = [1, 0, 2]

        result = tilde_U(U, kappa)

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

        result = tilde_U(U, kappa)

        @test result[1, :] == U[2, :]
        @test result[2, :] == U[1, :]
    end

    @testset "Edge cases" begin
        # Empty matrix
        U = Matrix{Float64}(undef, 0, 0)
        kappa = Int[]
        result = tilde_U(U, kappa)
        @test size(result) == (0, 0)

        # Single element
        U = reshape([1.0], 1, 1)
        kappa = [1]
        result = tilde_U(U, kappa)
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

        result = tilde_U(U, kappa)

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
        @test_throws BoundsError tilde_U(U, [3, 1])  # Index 3 is out of bounds

        # Test mismatched dimensions
        @test_throws DimensionMismatch tilde_U(U, [1, 2, 3])
    end

    @testset "Type stability" begin
        U = [
            1.0 2.0
            3.0 4.0
        ]
        kappa = [1, 2]

        # Test that output type matches input type
        result = tilde_U(U, kappa)
        @test eltype(result) == eltype(U)

        # Test with different types
        U_int = [1 2; 3 4]
        result_int = tilde_U(U_int, kappa)
        @test eltype(result_int) == eltype(U_int)
    end
end

@testset "Z function tests" begin
    # Test Néel state configuration
    # Sites: 0-1-2 in a line, each connected to neighbors
    nn = [(1, 2), (2, 3), (1, 3)]  # nearest neighbor bonds
    kappa_up = [0, 1, 0]   # up spins at sites 1
    kappa_down = [1, 0, 2] # down spins at sites 1 and 3 

    @test KagomeDSL.Z(nn, kappa_up, kappa_down) == 2
end

@testset "reevaluateW! tests" begin
    @testset "Basic functionality" begin
        # Create a simple 2x2 test case
        U_up = [1.0 0.2; 0.2 1.0]
        U_down = [1.0 0.3; 0.3 1.0]
        ham = Hamiltonian{2,2}(
            U_up,
            U_down,
            zeros(4, 4),
            [],
            zeros(4, 3),
            zeros(4, 1),
        )
        kappa_up = [1, 2]
        kappa_down = [2, 1]
        mc = MC(ham, kappa_up, kappa_down, zeros(2, 2), zeros(2, 2))

        KagomeDSL.reevaluateW!(mc)

        tilde_U_up = tilde_U(U_up, kappa_up)
        tilde_U_down = tilde_U(U_down, kappa_down)
        U_upinvs = tilde_U_up \ I
        U_downinvs = tilde_U_down \ I
        W_up_expected = U_up * U_upinvs
        W_down_expected = U_down * U_downinvs

        @test isapprox(mc.W_up, W_up_expected, atol = 1e-10)
        @test isapprox(mc.W_down, W_down_expected, atol = 1e-10)
    end

    @testset "Matrix properties" begin
        # Test with larger matrices
        n = 4
        U_up = Matrix{Float64}(I, n, n) + 0.1 * rand(n, n)
        U_down = Matrix{Float64}(I, n, n) + 0.1 * rand(n, n)
        U_up = (U_up + U_up') / 2  # Make symmetric
        U_down = (U_down + U_down') / 2

        ham = Hamiltonian{4,4}(
            U_up,
            U_down,
            zeros(16, 16),
            [],
            zeros(16, 5),
            zeros(16, 3),
        )
        kappa_up = [1, 2, 3, 4]
        kappa_down = [4, 3, 2, 1]
        mc = MC(ham, kappa_up, kappa_down, zeros(n, n), zeros(n, n))

        KagomeDSL.reevaluateW!(mc)

        # Test matrix dimensions
        @test size(mc.W_up) == (n, n)
        @test size(mc.W_down) == (n, n)

        # Test symmetry preservation
        @test isapprox(mc.W_up, mc.W_up', atol = 1e-10)
        @test isapprox(mc.W_down, mc.W_down', atol = 1e-10)
    end

end

@testset "update_W! and update_W_matrices tests" begin
    @testset "update_W! basic functionality" begin
        # Simple 2x2 test case
        W = [1.0 0.2; 0.2 1.0]
        W_original = copy(W)

        # Test update with l=1, K=2
        col_cache = Vector{eltype(W)}(undef, size(W, 1))
        row_cache = Vector{eltype(W)}(undef, size(W, 2))
        update_W!(W, 1, 2, col_cache, row_cache)

        # Verify elements manually
        factor = W[1, 1] / W[2, 1]
        @test isapprox(W[1, 1], W_original[1, 1] - factor * (W_original[2, 1] - 1.0))
        @test isapprox(W[1, 2], W_original[1, 2] - factor * W_original[2, 2])

        # Check original matrix wasn't modified
        @test W != W_original
    end

    @testset "update_W! matrix properties" begin
        n = 4
        W = Matrix{ComplexF64}(I, n, n) + 0.1 * rand(n, n)
        W = (W + W') / 2  # Make symmetric
        W_original = copy(W)

        # Test with random l and K
        l = rand(1:n)
        K = rand(1:n)
        while K == l  # Ensure K ≠ l
            K = rand(1:n)
        end
        col_cache = Vector{ComplexF64}(undef, size(W, 1))
        row_cache = Vector{ComplexF64}(undef, size(W, 2))
        update_W!(W, l, K, col_cache, row_cache)

        # Test matrix dimensions preserved
        @test size(W) == size(W_original)

        # Test no NaN or Inf values
        @test !any(isnan, W)
        @test !any(isinf, W)
    end

    @testset "update_W_matrices" begin
        # Create test MC object
        n = 3
        W_up = Matrix{Float64}(I, n, n) + 0.1 * rand(n, n)
        W_down = Matrix{Float64}(I, n, n) + 0.1 * rand(n, n)
        W_up = (W_up + W_up') / 2
        W_down = (W_down + W_down') / 2

        # Create mock MC struct with necessary fields
        mc = MC(
            Hamiltonian{3,3}(
                zeros(n, n),
                zeros(n, n),
                zeros(n^2, n^2),
                [],
                zeros(n, 4),
                zeros(n, 2),
            ),
            zeros(Int, n),
            zeros(Int, n),
            copy(W_up),
            copy(W_down),
        )

        # Store original matrices
        W_up_original = copy(mc.W_up)
        W_down_original = copy(mc.W_down)

        # Test update
        l_up, K_up = 1, 2
        l_down, K_down = 2, 3

        update_W_matrices!(mc; K_up = K_up, K_down = K_down, l_up = l_up, l_down = l_down)

        # Verify both matrices were updated
        @test !isapprox(mc.W_up, W_up_original)
        @test !isapprox(mc.W_down, W_down_original)

        # Verify no NaN or Inf values
        @test !any(isnan, mc.W_up)
        @test !any(isnan, mc.W_down)
        @test !any(isinf, mc.W_up)
        @test !any(isinf, mc.W_down)
    end

    @testset "update_W! numerical stability" begin
        # Test with near-singular matrix
        W = [1.0 1e-10; 1e-10 1.0]

        # Should handle small values without producing NaN/Inf
        col_cache = Vector{eltype(W)}(undef, size(W, 1))
        row_cache = Vector{eltype(W)}(undef, size(W, 2))
        update_W!(W, 1, 2, col_cache, row_cache)
        @test !any(isnan, W)
        @test !any(isinf, W)

        # Test with larger matrix containing small values
        n = 5
        W = Matrix{Float64}(I, n, n) + 1e-10 * rand(n, n)
        rng = Random.Xoshiro(42)  # Fixed seed for reproducibility
        l = rand(rng, 1:n)
        K = rand(rng, 1:n)
        while K == l
            K = rand(rng, 1:n)
        end
        col_cache = Vector{eltype(W)}(undef, size(W, 1))
        row_cache = Vector{eltype(W)}(undef, size(W, 2))
        update_W!(W, l, K, col_cache, row_cache)
        @test !any(isnan, W)
        @test !any(isinf, W)
    end
end

@testset "is_occupied" begin
    kappa = [1, 0, 2, 0]
    @test is_occupied(kappa, 1) == true
    @test is_occupied(kappa, 2) == false
    @test is_occupied(kappa, 3) == true
    @test is_occupied(kappa, 4) == false
    @test_throws BoundsError is_occupied(kappa, 5)
end


@testset "update_configurations!" begin
    # Setup initial MC state
    n = 4
    N_up = 2
    N_down = 2
    W_up = Matrix{ComplexF64}(undef, n, N_up)
    W_down = Matrix{ComplexF64}(undef, n, N_down)
    W_up .= 1.0 + 0.0im
    W_down .= 1.0 + 0.0im
    mc = MC(
        Hamiltonian(N_up, N_down, DoubleKagome(1.0, n, n, (false, false))),
        [1, 0, 2, 0],
        [0, 1, 0, 2],
        copy(W_up),
        copy(W_down),
        Vector{ComplexF64}(undef, n),
        Vector{ComplexF64}(undef, N_up),
        Vector{ComplexF64}(undef, n),
        Vector{ComplexF64}(undef, N_down),
    )

    @testset "flag = 1" begin
        mc_copy = deepcopy(mc)
        update_configurations!(mc_copy, 1, 2, 1, 1, 1)

        # Test kappa_up update
        @test mc_copy.kappa_up[1] == 1
        @test mc_copy.kappa_up[2] == 0

        # Test kappa_down update
        @test mc_copy.kappa_down[1] == 0
        @test mc_copy.kappa_down[2] == 1
    end

    @testset "flag = 2" begin
        mc_copy = deepcopy(mc)
        update_configurations!(mc_copy, 2, 1, 2, 1, 1)

        # Test kappa_up update
        @test mc_copy.kappa_up[1] == 1
        @test mc_copy.kappa_up[2] == 0
        # Test kappa_down update
        @test mc_copy.kappa_down[1] == 0
        @test mc_copy.kappa_down[2] == 1
    end
end

@testset "Carlo.init!" begin
    param = Dict(:n1 => 2, :n2 => 2, :PBC => (false, false), :N_up => 6, :N_down => 6)
    mc = MC(param)
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 3, :seed => 123, :thermalization => 10),
    )

    # Test successful initialization with QR method
    Carlo.init!(mc, ctx, param)

    # Verify that tilde_U matrices are non-singular after initialization
    tilde_U_up = tilde_U(mc.Ham.U_up, mc.kappa_up)
    @test abs(det(tilde_U_up)) != 0.0

    N_down = param[:N_down]
    if N_down > 0
        tilde_U_down = tilde_U(mc.Ham.U_down, mc.kappa_down)
        @test abs(det(tilde_U_down)) != 0.0
    end
end

@testset "Carlo.init! singular configuration handling" begin
    n1 = 2
    n2 = 2
    ns = n1 * n2 * 3 # 12 sites
    N_up = 6
    N_down = 6

    # Create a mock U_up matrix that will lead to a singular tilde_U_up
    # Make all rows identical to guarantee singularity of tilde_U_up
    U_up_singular = rand(1, N_up) # Create one random row
    U_up_singular = repeat(U_up_singular, ns, 1) # Repeat it for all ns rows

    # Create a mock Hamiltonian with the singular U_up
    mock_ham_singular = Hamiltonian{N_up,N_down}(
        U_up_singular,
        rand(ns, N_down),
        zeros(ComplexF64, ns, ns),
        Tuple{Int,Int}[],
        zeros(ComplexF64, ns, N_up + 1),
        zeros(ComplexF64, ns, N_down - 1),
    )

    mc_singular = MC(
        mock_ham_singular,
        zeros(Int, ns),
        zeros(Int, ns),
        zeros(ns, N_up),
        zeros(ns, N_down),
    )

    ctx_singular = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 3, :seed => 123, :thermalization => 10),
    )

    # Expect an error to be thrown
    @test_throws ErrorException(
        "QR-based configuration is singular. The Hamiltonian may be rank-deficient.",
    ) Carlo.init!(
        mc_singular,
        ctx_singular,
        Dict(
            :n1 => n1,
            :n2 => n2,
            :PBC => (false, false),
            :N_up => N_up,
            :N_down => N_down,
        ),
    )
end

@testset "Carlo.sweep! and Carlo.measure!" begin
    param = Dict(:n1 => 2, :n2 => 2, :PBC => (false, false), :N_up => 6, :N_down => 6)
    mc = MC(param)

    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 3, :seed => 123, :thermalization => 10),
    )
    Carlo.init!(mc, ctx, param)

    # Test that a sweep runs without error
    @test begin
        Carlo.sweep!(mc, ctx)
        true
    end

    # Test that a measurement runs without error
    ctx.sweeps = 1
    @test begin
        Carlo.measure!(mc, ctx)
        true
    end
end

@testset "Checkpointing" begin
    param = Dict(:n1 => 2, :n2 => 2, :PBC => (false, false), :N_up => 3, :N_down => 3)
    mc1 = MC(param)
    mc2 = MC(param)

    # Create a temporary file for checkpointing
    mktemp() do path, io
        close(io) # Close the handle to allow HDF5 to use the file
        h5open(path, "w") do file
            group = create_group(file, "test")
            Carlo.write_checkpoint(mc1, group)
        end

        h5open(path, "r") do file
            group = file["test"]
            Carlo.read_checkpoint!(mc2, group)
        end
    end


    @test mc1.kappa_up == mc2.kappa_up
    @test mc1.kappa_down == mc2.kappa_down
end

using BitBasis

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
    x = LongBitStr(vcat(fill(1, 1), fill(0, 35), fill(0, 1), fill(1, 35)))
    # Sx Sy will flip 1 to spin down, the nearest neighbor will be spin up
    # this gives 2 configurations
    xprime = KagomeDSL.getxprime(ham, x)
    @test length(keys(xprime)) == 3
    # Sz interaction
    Sz_sum = (length(DK.nn) - 2) * (1 / 2)
    Sz_sum += 2 * (2 * (-1 / 4))
    @test xprime[(-1, -1, -1, -1)] == Sz_sum
    # Sx Sy interaction, only happens in the bonds having site 1
    Sx_Sy_sum = 0.0
    # Sx Sy = 1/2()
    # x1 is flip 1 down, filp 2 up
    x1 = LongBitStr(vcat([0], [1], fill(0, 34), [1], [0], fill(1, 34)))
    k_up = 2
    l_up = 1
    k_down = 1
    l_down = 1
    conf = (k_up, l_up, k_down, l_down)
    @test conf in keys(xprime)
    @test xprime[conf] == 1 / 2 * 2
    # x2 is flip 1 up, flip 3 down
    x2 = LongBitStr(vcat([0], [0], [1], fill(0, 33), [1], [1], [0], fill(1, 33)))
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
    up = BitVector(vcat([1], [0], fill(1, 5), fill(0, 5)))
    mc = KagomeDSL.MC(
        Dict(:n1 => 2, :n2 => 2, :PBC => (false, false), :N_up => 6, :N_down => 6),
    )
    # consider up: [1,0,1,1,1,1,1,0,0,0,0,0] down: [0,1,0,0,0,0,0,1,1,1,1,1]
    down = BitVector(vcat([0], [1], fill(0, 5), fill(1, 5)))
    @test KagomeDSL.getOL(mc, up, down) != 0.0
end

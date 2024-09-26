using BitBasis

@testset "Hamiltonian" begin
    DK = DoubleKagome(1.0, 4, 3, (false, false))
    H = KagomeDSL.Hmat(DK)
    # check H is Hermitian
    @test H ≈ H'
    @test isempty(findall(x -> !(x ≈ 0), H - H'))
    # test that the first row of H has only 2,3 position elements =1 and others =0
    @test H[1, 1] == 0
    @test H[1, 2] == 1
    @test H[1, 3] == 1
    for i = 4:36
        @test H[1, i] == 0
    end
    @test H[3, 13] == -1

    DK2 = KagomeDSL.DoubleKagome(1.0, 4, 3, (true, false))
    H2 = KagomeDSL.Hmat(DK2)
    @test isempty(findall(x -> !(x ≈ 0), H2 - H2'))
    @test H2 ≈ H2'

    @test H2[1, 2] == 1
    @test H2[1, 3] == 1
    @test CartesianIndex(1, 11) in DK2.nn
    @test H2[1, 11] == 1 # horizontal PBC
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
    x = LongBitStr(vcat(fill(1, 1), fill(0, 71)))
    xprime = KagomeDSL.getxprime(ham, x)
    @test length(keys(xprime)) == 1
    @test xprime[x] == 0.0

    x1 = LongBitStr(vcat([1], fill(0, 36), [1], fill(0, 34)))
    xprime1 = KagomeDSL.getxprime(ham, x1)
    @test length(keys(xprime1)) == 2
    @test xprime1[x1] == -0.5
    k1 = LongBitStr(vcat([0], [1], fill(0, 34), [1], fill(0, 35)))
    @test k1 in keys(xprime1)
    @test xprime1[k1] == 1.0 / 2.0 * 2
end

@testset "getOL" begin
    DK = DoubleKagome(1.0, 2, 2, (false, false))
    ham = KagomeDSL.Hamiltonian(6, 6, DK)
    # consider up: [1,0,1,1,1,1,1,0,0,0,0,0] down: [1,0,1,1,1,1,1,0,0,0,0,0]
    up = BitVector(vcat([1], [0], fill(1, 5), fill(0, 5)))
    @test KagomeDSL.getOL(ham, up, up) == 0.0
    # consider up: [1,0,1,1,1,1,1,0,0,0,0,0] down: [0,1,0,0,0,0,0,1,1,1,1,1]
    down = BitVector(vcat([0], [1], fill(0, 5), fill(1, 5)))
    @test KagomeDSL.getOL(ham, up, down) != 0.0
end


@testset "Gutzwiller" begin
    x = LongBitStr(vcat(fill(1, 36), fill(0, 36)))
    @test KagomeDSL.Gutzwiller(x) == 1.0
    x = LongBitStr(vcat(fill(1, 36), fill(0, 35), [1]))
    @test KagomeDSL.Gutzwiller(x) == 0.0
end

using KagomeDSL
using Test

@testset "Hamiltonian" begin
    DK = DoubleKagome(1.0, 4, 3, (false, false))
    H = KagomeDSL.Hmat(DK, 1.0, 18, 18)
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

    DK2 = DoubleKagome(1.0, 4, 3, (true, false))
    H2 = KagomeDSL.Hmat(DK2, 1.0, 18, 18)
    @test isempty(findall(x -> !(x ≈ 0), H2 - H2'))
    @test H2 ≈ H2'

    @test H2[1, 2] == 1
    @test H2[1, 3] == 1
    @test CartesianIndex(1, 11) in nearestNeighbor(DK2)
    @test H2[1, 11] == 1 # horizontal PBC
    @test H2[1, 3*4+1] == 0 # vertical OBC
end

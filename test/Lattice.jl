using KagomeDSL
using Test

@testset "Lattice" begin
    # 2,2 is a corner case, don't use it
    K1 = Kagome(1.0, 3, 3, (false, false))
    @test size(K1.distance_matrix) == (27, 27)
    # should always read i <= j for distance_matrix[i,j]

    # test the diagonal elements are all zeros
    for i = 1:27
        @test K1.distance_matrix[i, i] == 0.0
    end

    @test K1.distance_matrix[1, 2] ≈ 1.0
    @test K1.distance_matrix[1, 3] ≈ 1.0
    @test K1.distance_matrix[1, 4] ≈ 2.0
    @test K1.distance_matrix[1, 5] ≈ 3.0
    @test K1.distance_matrix[1, 6] ≈ √(2.5^2 + (√3 / 2)^2)
    @test K1.distance_matrix[1, 7] ≈ 4.0
    @test K1.distance_matrix[1, 10] ≈ 2.0
    @test K1.distance_matrix[1, 13] ≈ 2√3

    K2 = Kagome(1.0, 3, 3, (true, false))
    @test size(K2.distance_matrix) == (27, 27)

    # test the diagonal elements are all zeros
    for i = 1:27
        @test K2.distance_matrix[i, i] == 0.0
    end

    @test K2.distance_matrix[1, 2] ≈ 1.0
    @test K2.distance_matrix[1, 3] ≈ 1.0
    @test K2.distance_matrix[1, 4] ≈ 2.0
    @test K2.distance_matrix[1, 5] ≈ 3.0
    @test K2.distance_matrix[1, 6] ≈ √(2.5^2 + (√3 / 2)^2)
    @test K2.distance_matrix[1, 7] ≈ 2.0 # periodic boundary condition
    @test K2.distance_matrix[1, 8] ≈ 1.0
    @test K2.distance_matrix[1, 9] ≈ √3

end

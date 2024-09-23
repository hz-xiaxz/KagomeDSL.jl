using KagomeDSL
using Test

@testset "MC" begin
    param = Dict(
        :n1 => 4,
        :n2 => 3,
        :PBC => (true, false),
        :χ => 1.0,
        :N_up => 18,
        :N_down => 18,
    )
    test_mc = MC(param)
    for i = 1:36
        @test test_mc.conf_up[i] ⊻ test_mc.conf_down[i] == true
    end
end

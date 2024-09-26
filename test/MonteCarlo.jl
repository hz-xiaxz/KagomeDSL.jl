using KagomeDSL
using Test

@testset "MC" begin
    param = Dict(
        :n1 => 4,
        :n2 => 3,
        :PBC => (true, false),
        :N_up => 18,
        :N_down => 18,
    )
    test_mc = MC(param)
    @test sum(test_mc.conf_up) == 18
    @test sum(test_mc.conf_down) == 18
end

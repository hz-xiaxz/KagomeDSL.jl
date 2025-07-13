using KagomeDSL
using Test
using Random
using Carlo

@testset "S_xy correlator" begin
    params = Dict(
        :n1 => 2,
        :n2 => 2,
        :PBC => (true, true),
        :N_up => 6,
        :N_down => 6,
    )
    mc = MC(params)
    ns = params[:n1] * params[:n2] * 3
    ctx = Carlo.MCContext{Random.Xoshiro}(
        Dict(:binsize => 3, :seed => 123, :thermalization => 10),
    )
    Carlo.init!(mc, ctx, params)

    # Find a pair of sites (i, j) with opposite spins
    i_up, j_down = -1, -1
    for i in 1:ns
        if mc.kappa_up[i] != 0
            for j in 1:ns
                if mc.kappa_down[j] != 0
                    i_up = i
                    j_down = j
                    break
                end
            end
        end
        if i_up != -1
            break
        end
    end

    @test i_up != -1 && j_down != -1

    # --- Test the S-_i S+_j term ---
    # Manual calculation
    l_up = mc.kappa_up[i_up]
    l_down = mc.kappa_down[j_down]
    manual_corr = 0.5 * mc.W_up[j_down, l_up] * mc.W_down[i_up, l_down]

    # Calculation using spinInteraction!
    xprime = Dict{KagomeDSL.ConfigKey,Float64}()
    spinInteraction!(xprime, mc.kappa_up, mc.kappa_down, i_up, j_down)
    
    refactored_corr = 0.0
    for (conf, coeff) in pairs(xprime)
        # We are only interested in the S-_i S+_j term for this comparison
        # This corresponds to the configuration where up spin moves from i->j
        # and down spin moves from j->i.
        # In spinInteraction!, this is the second case: K_up=j, K_down=i
        if conf.K_up == j_down && conf.K_down == i_up
            update_up = mc.W_up[conf.K_up, conf.l_up]
            update_down = mc.W_down[conf.K_down, conf.l_down]
            ratio = update_up * update_down
            refactored_corr += -coeff * ratio
        end
    end
    
    @test manual_corr ≈ refactored_corr
end

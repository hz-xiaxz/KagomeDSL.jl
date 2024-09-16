abstract type AbstractAnsatz end

"""
    fast_update(U::AbstractMatrix, Uinvs::AbstractMatrix, newconf::BitStr{N,T}, oldconf::BitStr{N,T}) where {N,T}

Fast computing technique from Becca and Sorella 2017
"""
@inline function fast_update(
    U::AbstractMatrix,
    Uinvs::AbstractMatrix,
    newconf::Union{SubDitStr{D,N1,T1},DitStr{D,N1,T1}},
    oldconf::BitStr{N,T},
) where {D,N1,N,T,T1}
    @assert length(newconf) == N "The length of the new configuration should be the same as the old configuration, got: $(length(newconf))(old) and $N(new)"
    Rl = -1 # if not found should return error
    k = -1
    flag = 0
    @inbounds for i = 1:N
        if getindex(oldconf, i) == 1 && getindex(newconf, i) == 0
            Rl = i # the old position of the l-th electron
            flag += 1
        end
        if getindex(newconf, i) == 1 && getindex(oldconf, i) == 0
            k = i # the new position of the l-th electron, K = R_l'
            flag += 1
        end
        if flag == 2
            break
        end
    end
    if flag == 0
        return 1.0
    end
    l = sum(oldconf[1:Rl]) # l-th electron
    ratio = sum(U[k, :] .* Uinvs[:, l])
    return ratio
end

@doc raw"""

The observable ``O_L = \frac{<x|H|\psi_G>}{<x|\psi_G>}``
"""
@inline function getOL(
    U_up::AbstractMatrix,
    U_down::AbstractMatrix,
    conf_up::BitVector,
    conf_down::BitVector,
)
    conf = LongBitStr(vcat(conf_up, conf_down))
    L = length(conf) ÷ 2
    OL = 0.0
    U_upinvs = U_up[conf_up, :] \ I # do invs more efficiently
    U_downinvs = U_down[conf_down, :] \ I
    xprime = getxprime(orb, conf)
    conf_downstr = LongBitStr(conf_down)
    conf_upstr = LongBitStr(conf_up)
    @inbounds for (confstr, coff) in pairs(xprime)
        if confstr == conf
            OL += coff
        else
            OL +=
                coff *
                fast_update(U_up, U_upinvs, SubDitStr(confstr, 1, L), conf_upstr) *
                fast_update(U_down, U_downinvs, SubDitStr(confstr, L + 1, 2L), conf_downstr)
        end
    end
    return OL
end

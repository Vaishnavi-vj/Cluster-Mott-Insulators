ind(i, orbname) = (up=i + (2 * (dict[orbname] - 1)), down=i + (2 * (dict[orbname] - 1)) + 1)
neigh(i, norbs) = i + (2 * norbs)

dict_lop = Dict(L_ind => m_ind for (L_ind, m_ind) in zip(1:3, ['b', 'c', 'a']))
dict_jop = Dict(L_type => L_sign for (L_type, L_sign) in zip(["real", "eff"], [1, -1]))

########################################################################3
#
#                  SPIN ORBIT COUPLING
#
########################################################################

function socrot(N, norbs, state_list; ltype::String="eff")

    Lsign = dict_jop[ltype]

    M = length(state_list)
    H = spzeros(Complex{Float64}, M, M)

    """
    x
    if (abit[ind(i, 'b').up] != abit[ind(i, 'c').down])
    if (abit[ind(i, 'b').down] != abit[ind(i, 'c').up])

    y
    if (abit[ind(i, 'c').up] != abit[ind(i, 'a').down])
    if (abit[ind(i, 'c').down] != abit[ind(i, 'a').up])

    z
    if (abit[ind(i, 'a').up] != abit[ind(i, 'b').up])
    if (abit[ind(i, 'a').down] != abit[ind(i, 'b').down])
    """
    for a in 1:M
        abit = digits(state_list[a], base=2, pad=2 * N * norbs)
        for i in 1:2*norbs:2*norbs*N
            # SOC_x
            spinPairsX = [(ind(i, 'b').up, ind(i, 'c').down)]
            orbitalTerms!(H, i, i, abit, -im * (-1)^abit[ind(i, 'b').up] * Lsign, a, state_list)

            spinPairsX = [(ind(i, 'b').down, ind(i, 'c').up)]
            orbitalTerms!(H, i, i, abit, -im * (-1)^abit[ind(i, 'b').down] * Lsign, a, state_list)

            """
            if (abit[ind(i, 'b').up] != abit[ind(i, 'c').down])

                s_star = flip(abit, ind(i, 'b').up, ind(i, 'c').down)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'b').up+1:ind(i, 'c').down-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'b').up]) * jw_str * Lsign
                end
            end
            if (abit[ind(i, 'b').down] != abit[ind(i, 'c').up])

                s_star = flip(abit, ind(i, 'b').down, ind(i, 'c').up)
                jw_str = 1 #NOTE:WHY??
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'b').down]) * jw_str * Lsign
                end
            end
            """

            #SOC_y
            if (abit[ind(i, 'c').up] != abit[ind(i, 'a').down])

                s_star = flip(abit, ind(i, 'c').up, ind(i, 'a').down)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').down+1:ind(i, 'c').up-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] - 0.5 * jw_str * Lsign
                end
            end

            if (abit[ind(i, 'c').down] != abit[ind(i, 'a').up])

                s_star = flip(abit, ind(i, 'c').down, ind(i, 'a').up)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'c').down-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * jw_str * Lsign
                end
            end

            #SOC_z                   
            if (abit[ind(i, 'a').up] != abit[ind(i, 'b').up])

                s_star = flip(abit, ind(i, 'a').up, ind(i, 'b').up)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'b').up-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'a').up]) * jw_str * Lsign
                end
            end

            if (abit[ind(i, 'a').down] != abit[ind(i, 'b').down])

                s_star = flip(abit, ind(i, 'a').down, ind(i, 'b').down)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').down+1:ind(i, 'b').down-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'b').down]) * jw_str * Lsign
                end
            end
        end
    end
    return H
end

###################################################################################


function socrot(N, norbs, state_list, rotmats; ltype::String="eff")

    M = length(state_list)
    H = spzeros(Complex{Float64}, M, M)
    Lsign = dict_jop[ltype]

    for a in 1:M

        abit = digits(state_list[a], base=2, pad=2 * N * norbs)
        site = 1

        for i in 1:2*norbs:2*norbs*N

            # SOC_LzSx
            if (abit[ind(i, 'a').up] != abit[ind(i, 'b').down])

                s_star = flip(abit, ind(i, 'a').up, ind(i, 'b').down)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'b').down-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * rotmats[site][1, 1] * ((-1)^abit[ind(i, 'a').up]) * jw_str * Lsign
                end
            end

            if (abit[ind(i, 'a').down] != abit[ind(i, 'b').up])

                s_star = flip(abit, ind(i, 'a').down, ind(i, 'b').up)
                jw_str = 1#(-1)^sum(abit[p] for p in ind(i,'a').down+1:ind(i,'b').up-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * rotmats[site][1, 1] * ((-1)^abit[ind(i, 'a').down]) * jw_str * Lsign
                end

            end

            # SOC_LzSy
            if (abit[ind(i, 'a').up] != abit[ind(i, 'b').down])

                s_star = flip(abit, ind(i, 'a').up, ind(i, 'b').down)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'b').down-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * rotmats[site][1, 2] * jw_str * Lsign
                end
            end

            if (abit[ind(i, 'b').up] != abit[ind(i, 'a').down])

                s_star = flip(abit, ind(i, 'b').up, ind(i, 'a').down)
                jw_str = 1#(-1)^sum(abit[p] for p in ind(i,'a').down+1:ind(i,'c').down-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] - 0.5 * rotmats[site][1, 2] * jw_str * Lsign
                end

            end

            #SOC_LzSz                   
            if (abit[ind(i, 'a').up] != abit[ind(i, 'b').up])

                s_star = flip(abit, ind(i, 'a').up, ind(i, 'b').up)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'b').up-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * rotmats[site][1, 3] * ((-1)^abit[ind(i, 'a').up]) * jw_str * Lsign
                end
            end

            if (abit[ind(i, 'a').down] != abit[ind(i, 'b').down])

                s_star = flip(abit, ind(i, 'a').down, ind(i, 'b').down)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').down+1:ind(i, 'b').down-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * rotmats[site][1, 3] * ((-1)^abit[ind(i, 'b').down]) * jw_str * Lsign
                end
            end


            # SOC_LySx
            if (abit[ind(i, 'c').up] != abit[ind(i, 'a').down])

                s_star = flip(abit, ind(i, 'c').up, ind(i, 'a').down)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').down+1:ind(i, 'c').up-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * rotmats[site][2, 1] * ((-1)^abit[ind(i, 'c').up]) * jw_str * Lsign
                end
            end

            if (abit[ind(i, 'c').down] != abit[ind(i, 'a').up])

                s_star = flip(abit, ind(i, 'a').up, ind(i, 'c').down)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'c').down-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * rotmats[site][2, 1] * ((-1)^abit[ind(i, 'c').down]) * jw_str * Lsign
                end
            end

            #SOC_LySy
            if (abit[ind(i, 'c').up] != abit[ind(i, 'a').down])

                s_star = flip(abit, ind(i, 'c').up, ind(i, 'a').down)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').down+1:ind(i, 'c').up-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * rotmats[site][2, 2] * jw_str * Lsign
                end
            end

            if (abit[ind(i, 'c').down] != abit[ind(i, 'a').up])

                s_star = flip(abit, ind(i, 'c').down, ind(i, 'a').up)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'c').down-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] - 0.5 * rotmats[site][2, 2] * jw_str * Lsign
                end
            end

            # SOC_LySz
            if (abit[ind(i, 'c').up] != abit[ind(i, 'a').up])

                s_star = flip(abit, ind(i, 'c').up, ind(i, 'a').up)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'c').up-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * rotmats[site][2, 3] * ((-1)^abit[ind(i, 'c').up]) * jw_str * Lsign
                end
            end

            if (abit[ind(i, 'c').down] != abit[ind(i, 'a').down])

                s_star = flip(abit, ind(i, 'a').down, ind(i, 'c').down)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').down+1:ind(i, 'c').down-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * rotmats[site][2, 3] * ((-1)^abit[ind(i, 'a').down]) * jw_str * Lsign
                end

            end

            # SOC_LxSx
            if (abit[ind(i, 'b').up] != abit[ind(i, 'c').down])

                s_star = flip(abit, ind(i, 'b').up, ind(i, 'c').down)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'b').up+1:ind(i, 'c').down-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * rotmats[site][3, 1] * ((-1)^abit[ind(i, 'b').up]) * jw_str * Lsign
                end
            end

            if (abit[ind(i, 'b').down] != abit[ind(i, 'c').up])

                s_star = flip(abit, ind(i, 'b').down, ind(i, 'c').up)
                jw_str = 1
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * rotmats[site][3, 1] * ((-1)^abit[ind(i, 'b').down]) * jw_str * Lsign
                end
            end


            # SOC_LxSy
            if (abit[ind(i, 'b').up] != abit[ind(i, 'c').down])

                s_star = flip(abit, ind(i, 'b').up, ind(i, 'c').down)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'b').up+1:ind(i, 'c').down-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * rotmats[site][3, 2] * jw_str * Lsign
                end
            end

            if (abit[ind(i, 'b').down] != abit[ind(i, 'c').up])

                s_star = flip(abit, ind(i, 'b').down, ind(i, 'c').up)
                jw_str = 1#(-1)^sum(abit[p] for p in ind(i,'b').down+1:ind(i,'c').up-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] - 0.5 * rotmats[site][3, 2] * jw_str * Lsign
                end
            end

            # SOC_LxSz
            if (abit[ind(i, 'b').up] != abit[ind(i, 'c').up])

                s_star = flip(abit, ind(i, 'b').up, ind(i, 'c').up)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'b').up+1:ind(i, 'c').up-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * rotmats[site][3, 3] * ((-1)^abit[ind(i, 'b').up]) * jw_str * Lsign
                end
            end

            if (abit[ind(i, 'b').down] != abit[ind(i, 'c').down])

                s_star = flip(abit, ind(i, 'b').down, ind(i, 'c').down)
                jw_str = (-1)^sum(abit[p] for p in ind(i, 'b').down+1:ind(i, 'c').down-1)
                b = findstate(s_star[2], M, state_list)

                if (b > 0)
                    H[a, b] = H[a, b] + 0.5 * im * rotmats[site][3, 3] * ((-1)^abit[ind(i, 'c').down]) * jw_str * Lsign
                end
            end

            site = site + 1

        end
    end

    return H

end

################################################################################################

function crystalFieldSplitting!(H, i, abit, a, delta, state_list, orbitals)
    for m in orbitals
        for n in orbitals
            (m == n) && continue
            interOrbitalTerms!(H, i, i, abit, m, n, delta, a, state_list)
        end
    end
end

function CF_splitting(N, norbs, delta, state_list)

    M = length(state_list)
    H = spzeros(Complex{Float64}, M, M)
    orbitals = 'a':'c'
    for a in 1:M
        abit = digits(state_list[a], base=2, pad=2 * N * norbs)
        for i in 1:2*norbs:2*norbs*N
            crystalFieldSplitting!(H, i, abit, a, delta, state_list, orbitals)
        end
    end
    return H
end


function CF_splitting_deltasite(N, norbs, delta, state_list)

    M = length(state_list)
    Hcf = spzeros(Complex{Float64}, M, M)
    orbitals = 'a':'a'+norbs-1
    for a in 1:M
        abit = digits(state_list[a], base=2, pad=2 * N * norbs)
        site = 1
        for i in 1:2*norbs:2*norbs*N
            crystalFieldSplitting!(H, i, abit, a, delta[site], state_list, orbitals)
            site = site + 1
        end
    end
    return H
end
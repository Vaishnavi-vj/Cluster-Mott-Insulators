"""
--------------------------------------------------------------------------------
Hamiltonian functions:                      
Construct the Hubbard-Kanamori Hamiltonian for multi-orbital/site case,
with conservation of various symmetries as indicated            
--------------------------------------------------------------------------------
"""

function UTerm(U, i, abit, orbital)
    up, down = ind(i, orbital)
    return U * abit[up] * abit[down]
end

"""
Symmetric upon m <-> n
"""
function UprimeTerm(Up, Jh, i, abit, m, n)
    upₘ, downₘ = ind(i, m)
    upₙ, downₙ = ind(i, n)
    UprimeJ = 0.5 * (Up - Jh)
    return UprimeJ * (abit[upₘ] * abit[upₙ] + abit[downₘ] * abit[downₙ]) + Up * abit[upₘ] * abit[downₙ]
end

function spinExchange(jordanWignerString, Jh)
    return -0.5 * jordanWignerString * Jh
end

function pairHop(jordanWignerString, Jh)
    return 0.5 * jordanWignerString * Jh
end

function computeJordanWignerString(bitIndexPairs::Vector{Pair{Int64,Int64}})
    return (-1)^mapfoldl((x) -> x[1] & (1 << (x[2] - 1)), +, bitIndexPairs)
end

function spinExchangePairHop(i, a, abit, m, n, state_list, Jh)
    upₘ, downₘ = ind(i, m)
    upₙ, downₙ = ind(i, n)
    aupₘ, adownₘ = abit[[ind(i, m)...]]
    aupₙ, adownₙ = abit[[ind(i, n)...]]
    pairIndex = 0
    if !(aupₘ == aupₙ || adownₘ == adownₙ)
        s_star = flip(state_list[a], upₘ, upₙ, downₘ, downₙ)
        pairIndex = findstate(s_star, state_list)
        s_star = flip(first(flip(abit, upₘ, upₙ)), downₘ, downₙ)
    end
    if pairIndex > 0
        jordanWignerString = (-1)^(aupₘ + first(s_star)[upₙ])
        """
        Why is this only computed using up indexes?
        """
        #jordanWignerString = computeJordanWignerString((state_list[a], s_star) .=> (upₘ, upₙ))
        if (aupₘ + adownₘ == aupₙ + adownₙ)
            return spinExchange(jordanWignerString, Jh), pairIndex
        else
            return pairHop(jordanWignerString, Jh), pairIndex
        end
    else
        return zeros(Int64, 2)
    end
end

function interactionTerms!(H, state_list, a, abit, i, norbs, U, Up, Jh)
    orbitals = 'a':'a'+norbs-1
    for m in orbitals
        for n in orbitals
            if (m == n)
                #U-term
                H[a, a] += UTerm(U, i, abit, m)
            else
                # (U'-J) and U'-terms
                H[a, a] += UprimeTerm(Up, Jh, i, abit, m, n)

                # J-terms: spin-exchange and pair hopping
                Jterm, pairIndex = spinExchangePairHop(i, a, abit, m, n, state_list, Jh)
                if !(pairIndex == 0)
                    H[a, pairIndex] += Jterm
                end
            end
        end
    end
    return nothing
end

function spinHopping(i, j, referenceSite, currentSite, a, abit, state_list)
    if i == j
        @goto exitEarly
    end

    s_star = flip(state_list[a], i, j)
    pairIndex = findstate(s_star, state_list)

    if pairIndex > 0
        sumRange = referenceSite < currentSite ? range(i + 1, j - 1) : range(j + 1, i - 1)
        jordanWignerString = (-1)^sum(abit[sumRange])
        return jordanWignerString, pairIndex
    else
        @label exitEarly
        return zeros(Int64, 2)
    end
end

function orbitalTerms!(H, referenceSite, currentSite, spinPairs, abit, t, a, state_list)
    for (spinᵢ, spinⱼ) in spinPairs
        jordanWignerString, pairIndex = spinHopping(spinᵢ, spinⱼ, referenceSite, currentSite, a, abit, state_list)
        if !(jordanWignerString == 0)
            H[a, pairIndex] -= t * jordanWignerString
        end
    end
end

function intraOrbitalTerms!(H, referenceSite, currentSite, abit, orbital, tm, a, state_list)
    spinPairs = zip(ind(referenceSite, orbital), ind(currentSite, orbital))
    orbitalTerms!(H, referenceSite, currentSite, spinPairs, abit, tm, a, state_list)
end

function interOrbitalTerms!(H, referenceSite, currentSite, abit, m, n, tmn, a, state_list)
    spinPairs = zip(ind(referenceSite, m), ind(currentSite, n))
    orbitalTerms!(H, referenceSite, currentSite, spinPairs, abit, tmn, a, state_list)
end

function hoppingTerms!(H, tarr, a, abit, i, norbs, N, state_list)
    orbitals = 'a':'a'+norbs-1
    for (tm, tmn) in tarr
        j = (i != 2 * norbs * N - (2 * norbs - 1)) ? (neigh(i, norbs)) : 1 #Simplify later
        for m in orbitals
            for n in orbitals
                if m == n
                    intraOrbitalTerms!(H, i, j, abit, m, tm, a, state_list)
                else
                    interOrbitalTerms!(H, i, j, abit, m, n, tmn, a, state_list)
                end
            end
        end
    end
    return nothing
end
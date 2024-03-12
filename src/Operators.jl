
"""
ToDo: Move these globals elsewhere
"""
ind(i, orbname) = (
    up=i + (2 * (dict[orbname] - 1)), down=i + (2 * (dict[orbname] - 1)) + 1)
neigh(i, norbs) = i + (2 * norbs)

dict_lop = Dict(L_ind => m_ind for (L_ind, m_ind) in zip(1:3, ['b', 'c', 'a']))
dict_jop = Dict(L_type => L_sign for (L_type, L_sign) in zip(["real", "eff"], [1, -1]))

"""
--------------------------------------------------------------------------------
Operators list:
S_op, S_opi, Stot, Stot_tr
sum_L_op, sum_L_optr, L_opi, L_optri
J_op, Jtot
--------------------------------------------------------------------------------
"""

"""
--------------------------------------------------------------------------------
S operators
--------------------------------------------------------------------------------
"""

"""
(a) Without translation
"""

function S_op_term!(S, i, j, currentIndex, abit, term)
    abit[i] == abit[j] && return nothing
    s_star = flip(abit, i, j)
    pairIndex = findstate(s_star, state_list)
    if (pairIndex > 0)
        S[currentIndex, pairIndex] += term
    end
end

function S_op_term!(Sᵢ, i, m, a, abit, ::Val{'x'})
    index = ind(i, m)
    S_op_term!(Sᵢ, index.up, index.down, a, abit, 0.5)
end
function S_op_term!(Sᵢ, i, m, a, abit, ::Val{'y'})
    index = ind(i, m)
    S_op_term!(Sᵢ, index.up, index.down, a, abit, 0.5 * im * (-1)^abit[index.down])
end
function S_op_term!(Sᵢ, i, m, a, abit, ::Val{'z'})
    Sᵢ[a, a] += 0.5 * (abit[ind(i, m).up] - abit[ind(i, m).down])
end
function S_op_term!(Sᵢ, i, m, a, abit, ::Val{"all"})
    Sxi, Syi, Szi = Sᵢ
    index = ind(i, m)
    S_op_term!(Sxi, index.up, index.down, a, abit, 0.5)
    S_op_term!(Syi, index.up, index.down, a, abit, 0.5 * im * (-1)^abit[index.down])
    Szi[a, a] += 0.5 * (abit[index.up] - abit[index.down])
end

function S_op(N, norbs, state_list, choice)
    sitei = 1:(2*norbs):(2*norbs*N)
    return S_opi(N, norbs, state_list, sitei, choice)
end

function S_opi(N, norbs, state_list, sitei, choice)
    M = length(state_list)
    if choice == "all"
        Sᵢ = (
            spzeros(ComplexF64, M, M),
            spzeros(ComplexF64, M, M),
            spzeros(ComplexF64, M, M),
        )
    else
        Sᵢ = spzeros(ComplexF64, M, M)
    end

    orbitals = 'a':('a'+norbs-1)
    for a in 1:M
        abit = digits(state_list[a], base=2, pad=2 * N * norbs)
        for i in sitei
            for m in orbitals
                S_op_term!(Sᵢ, i, m, a, abit, Val(choice))
            end
        end
    end
    return Si
end

"""
(3) STOT functions
Performs the operation <ψ₀|Stot|ψ₀> = s(s+1), where
|ψ₀> is the groundstate eigenvector obtained from diagonalising
Hamiltonians from (2). Note that the symmetries taken into account
while obtaining |ψ₀> should also be considered while computing Stot.
Note: Stot has a structure similar to Heisenberg Hamiltonian
"""

"""
(3a) Stot function to calculate <ψ₀|Stot|ψ₀> when |ψ₀> is obtained from Hubbard_nfsz_multi
"""

function StotDiagonal!(S², i, m, a, abit, J)
    up, down = ind(i, m)
    bits = abit[up:down]
    if isone(sum(bits))
        S²[a, a] += J * (3 / 4)
    end
end

function StotOffDiagonal!(S², i, j, m, n, a, abit, J)
    upim, downim = ind(i, m)
    upjn, downjn = ind(j, n)
    bitsim = abit[upim:downim]
    bitsjn = abit[upjn:downjn]

    if (bitsim == bitsjn) && isone(sum(bitsim))
        S²[a, a] += J * (1 / 4)
    elseif (bitsim == reverse(bitsjn)) && isone(sum(bitsim))
        S²[a, a] -= J * (1 / 4)
        s_star = flip(a, upim, upjn, downim, downjn)
        pairIndex = findstate(s_star, state_list)

        if (pairIndex > 0)
            S²[a, pairIndex] += J * (1 / 2)
        end
    end
end
function StotOffDiagonalTr!(S², i, j, m, n, a, abit, J, Nconst, trfac, representativeFunc)
    upim, downim = ind(i, m)
    upjn, downjn = ind(j, n)
    bitsim = abit[upim:downim]
    bitsjn = abit[upjn:downjn]

    if (bitsim == bitsjn) && isone(sum(bitsim))
        S²[a, a] += J * (1 / 4)
    elseif (bitsim == reverse(bitsjn)) && isone(sum(bitsim))
        S²[a, a] -= J * (1 / 4)
        s_star = flip(a, upim, upjn, downim, downjn)
        pairInfo = representativeFunc(s_star, N, norbs)
        pairIndex = findstate(pairInfo[1], state_list)
        if (pairIndex > 0)
            S²[a, b] += 0.5 * J * (Nconst[pairIndex] / Nconst[a]) * (trfac^pairInfo[2])
        end
    end
end

function Stot(N, norbs, ψ₀, state_list)
    J = 1.0 #Question: Remove this J parameter?
    M = length(state_list)
    S² = spzeros(Complex{Float64}, M, M)

    sites = 1:(2*norbs):(2*norbs*N-1)
    sitePairs = Iterators.product(sites, sites)
    orbitals = 'a':('a'+norbs-1)
    orbitalPairs = Iterators.product(orbitals, orbitals)

    """
    Question:
    Isnt this calculation symmetrical under the condition of 
    i <-> j? 

    In which case we can cut the number of loops in half?
    """

    for a in 1:M
        abit = digits(state_list[a], base=2, pad=2 * N * norbs)
        for (i, j) in sitePairs
            for (m, n) in orbitalPairs
                # Diagonal case: same site and same orbital
                if (i == j) && (m == n)
                    StotDiagonal!(S, i, m, a, abit, J)
                else
                    # Off-diagonal cases
                    StotOffDiagonal!(S, i, j, m, n, a, abit, J)
                end
            end
        end
    end
    # Calculate and return <ψ₀|S^2|ψ₀> = s(s+1)
    return (ψ₀' * S² * ψ₀)
end

"""
(3b) Function to calculate <ψ₀|Stot|ψ₀> when |ψ₀> is obtained from Hubbard_tr_multi: i.e., with translations
"""

function Stot_tr(N::Int64,
    nf::Int64,
    k::Int64,
    norbs::Int64,
    ψ₀,
    state_list::Tuple{Array{Int64,1},Array{Float64,1}})
    J = 1.0
    rep_states = state_list[1]
    Nconst = state_list[2]

    M = length(rep_states)
    S2 = Complex.(spzeros(M, M))
    trfac = exp(-im * 2 * pi * k / N)

    if iseven(nf)
        representativeFunc = representative_even
    else
        representativeFunc = representative_odd
    end

    for a in 1:M
        abit = digits(rep_states[a], base=2, pad=2 * N * norbs)
        for i in 1:(2*norbs):(2*norbs*N-1)
            for j in 1:(2*norbs):(2*norbs*N-1)
                for m in 'a':('a'+norbs-1)
                    for n in 'a':('a'+norbs-1)
                        # The diagonal case where it is the same site and orbital
                        if (i == j && m == n)
                            StotDiagonal!(S2, i, m, a, abit, J)
                        else
                            StotOffDiagonalTr!(S2, i, j, m, n, a, abit, J, Nconst, trfac, representativeFunc)
                        end
                    end
                end
            end
        end
    end
    # Calculate and return <ψ₀|S^2|ψ₀> = s(s+1)
    return (ψ₀' * S2 * ψ₀)
end

"""
--------------------------------------------------------------------------------
L Operators
--------------------------------------------------------------------------------
"""

function jordanWignerString(reference, current, i, j, abit)
    sumRange = reference < current ? range(i + 1, j - 1) : range(j + 1, i - 1)
    return (-1)^sum(abit[sumRange])
end

function L_op_term!(Li, i, j, referenceOrbital, currentOrbital, currentIndex, abit, term)
    abit[i] == abit[j] && return nothing
    s_star = flip(abit, i, j)
    pairIndex = findstate(s_star[2], M, state_list)
    if (pairIndex > 0)
        if abit[j] == 0
            term = -term
        end
        Li[currentIndex, pairIndex] += term * jordanWignerString(referenceOrbital, currentOrbital, i, j, abit)
    end
end

function L_op_term_tr!(Li, i, j, referenceOrbital, currentOrbital, currentIndex, abit, Lsign, Nconst, trfac, repFunc)
    abit[i] == abit[j] && nothing
    s_star = flip(abit, i, j)
    pairInfo = repFunc(s_star[2], s_star[1], N, norbs)
    pairIndex = findstate(pairInfo[1], M, rep_states)

    if (pairIndex > 0)
        if (abit[i] == 0)
            Lsign = -Lsign
        end
        im * Lsign * jw_str * b_info[3] * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
        Li[currentIndex, pairIndex] +=
            im * Lsign * pairInfo[3] * (Nconst[pairIndex] / Nconst[currentIndex]) * trfac^pairInfo[2] *
            jordanWignerString(referenceOrbital, currentOrbital, i, j, abit)
    end
end

"""
Question: Why are these called sum_L_op instead of L_op?
"""

function sum_L_op(N, norbs, state_list, choice; ltype::String="eff")
    sites = 1:(2*norbs):(2*norbs*N)
    return L_opi(N, norbs, state_list, choice, sites, ltype=ltype)
end

function sum_L_optr(N, norbs, nf, k, ψ₀, state_list, choice; ltype::String="eff")
    sites = 1:(2*norbs):(2*norbs*N)
    return L_optri(N, norb, nf, k, ψ₀, state_list, choice, sites, ltype=ltype)
end

function L_opi(N, norbs, state_list, choice, sitei; ltype::String="eff")
    M = length(state_list)
    Li = spzeros(ComplexF64, M, M)

    Lsign = dict_jop[ltype]
    m = dict_lop[choice]
    n = m != 'c' ? m + 1 : 'a' #Question: Why does this exist? -> To loop in the orbitals

    for a in 1:M
        abit = digits(state_list[a], base=2, pad=2 * N * norbs)
        for i in sitei
            L_op_term!(Li, ind(i, m).up, ind(i, n).up, m, n, a, abit, Lsign * im)
            L_op_term!(Li, ind(i, m).down, ind(i, n).down, m, n, a, abit, Lsign * im)
        end
    end
    return Li
end

function L_optri(N, norbs, nf, k, ψ₀, state_list, choice, sitei; ltype::String="eff")
    rep_states = state_list[1]
    Nconst = state_list[2]

    M = length(rep_states)
    Li = Complex.(spzeros(M, M))
    Lsign = dict_jop[ltype]
    trfac = exp((-im * 2 * pi * k) / N)

    m = dict_lop[choice]
    n = m != 'c' ? m + 1 : 'a'

    if iseven(nf)
        representativeFunc = representative_even
    else
        representativeFunc = representative_odd
    end

    for a in 1:M
        abit = digits(rep_states[a], base=2, pad=2 * N * norbs)
        for i in sitei
            L_op_term_tr!(Li, ind(i, m).up, ind(i, n).up, m, n, a, abit, Lsign, Nconst, trfac, representativeFunc)
            L_op_term_tr!(Li, ind(i, m).down, ind(i, n).down, m, n, a, abit, Lsign, Nconst, trfac, representativeFunc)
        end
    end
    return Li * ψ₀
end

"""
--------------------------------------------------------------------------------
J operators
--------------------------------------------------------------------------------
"""

function J_op(N, norbs, state_list; choice::String="eff")
    M = length(state_list)
    Lsign = dict_jop[choice]

    Sops = spzeros(ComplexF64, M, M, 3) #Possibly change 3->norbs
    Lops = spzeros(ComplexF64, M, M, 3) #Possibly change 3->norbs

    Sx = Sops[:, :, 1]
    Sy = Sops[:, :, 2]
    Sz = Sops[:, :, 3]

    Lx = Lops[:, :, 1]
    Ly = Lops[:, :, 2]
    Lz = Lops[:, :, 3]

    """
    Question: 
    Is is necessary we only have 3 Lop and 3 Sop operators? 
    """
    for (m, Lcurrent) in zip(('b', 'c', 'a'), (Lx, Ly, Lz))
        n = m != 'c' ? m + 1 : 'a'
        for a in 1:M
            abit = digits(state_list[a], base=2, pad=2 * N * norbs)
            for i in 1:(2*norbs):(2*norbs*N)
                #L terms
                L_op_term!(Lcurrent, ind(i, m).up, ind(i, n).up, m, n, a, abit, Lsign * im)
                L_op_term!(Lcurrent, ind(i, m).down, ind(i, n).down, m, n, a, abit, Lsign * im)

                #S terms
                S_op_term!((Sx, Sy, Sz), i, m, a, abit, Val("all"))
            end
        end
    end
    return mapreduce((x) -> (Sops[:, :, i] + Lops[:, :, i])^2, +, 1:3)
end

"""
Question: What does L2 and S2 mean?
"""
function Jtot(L2, S2, soc, ψ₀; ltype::String="eff")
    return ψ₀' * (L2 + S2 + 2 * dict_jop[ltype] * (soc)) * ψ₀
end

"""
--------------------------------------------------------------------------------
Density functions
--------------------------------------------------------------------------------
"""

function site_density(N, norbs, ψ, state_list, i)
    M = length(state_list)
    ρ = zeros(ComplexF64, M, M)
    orbitals = collect('a':('a'+norbs-1))
    for a in 1:M
        abit = digits(state_list[a], base=2, pad=2 * N * norbs)
        for m in eachindex(orbitals)
            up, down = ind(i, orbitals[m])
            ρi[a, a] += abit[up] + abit[down]
        end
    end
    return ρ * ψ
end

function orb_density(N, norbs, ψ, state_list)
    M = length(state_list)
    ρ = zeros(ComplexF64, M, M, norbs)
    orbitals = collect('a':('a'+norbs-1))
    sites = 1:(2*norbs):(2*norbs*N-1)
    for a in 1:M
        abit = digits(state_list[a], base=2, pad=2 * N * norbs)
        for i in sites, m in eachindex(orbitals)
            up, down = ind(i, orbitals[m])
            ρ[:, :, m] += abit[up] + abit[down]
        end
    end
    return mapreduce((x) -> ψ' * ρ[:, :, x] * ψ, +, 1:norbs, init=zeros(ComplexF64, M, M))
end
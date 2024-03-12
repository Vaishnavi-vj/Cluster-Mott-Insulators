include("ClusterHubbard_multiorbital.jl")

#"""
#Possible parameters type
#"""
#@kwdef mutable struct HamiltonianParameters
#    N::Int64
#    norbs::Int64
#    U::Float64 = zero(Float64)
#    U′::Float64 = zero(Float64)
#    Jₕ::Float64 = zero(Float64)
#    t::Vector{Tuple{Float64,Float64}} = Vector{Tuple{Float64,Float64}}(undef, 0)
#    states::Vector{Int64} = Vector{Int64}(undef, 0)
#    #symmetries::Vector{Function}
#    function HamiltonianParameters(N::Int64, norbs::Int64)
#        return new(; N=N, norbs=norbs)
#    end
#end

#"""
#Possible Hamiltonian generating function
#"""
#function createHamiltonain(parameters, functions)
#    state_list = parameters.states
#    M = length(state_list)
#    H = spzeros(ComplexF64, M, M)
#
#    for a in eachindex(state_list)
#        abit = digits(state_list[a], base=2, pad=2 * norbs * N)
#        for i in 1:2*norbs:2*norbs*N-1
#            for f in functions
#                f(parameters, abit, a, i)
#            end
#        end
#    end
#end

function Hubbard_nfsz_multi(N::Int64,
    norbs::Int64,
    U::Float64,
    Up::Float64,
    Jh::Float64,
    tarr::Array{Tuple{Float64,Float64},1},
    state_list::Array{Int64,1})

    M = length(state_list)
    H = spzeros(Complex, M, M)
    tarr = (N != 2) ? tarr : [0.5 .* i for i in tarr]

    for a in 1:M
        abit = digits(state_list[a], base=2, pad=2 * norbs * N)
        for i in 1:2*norbs:2*norbs*N-1
            interactionTerms!(H, state_list, a, abit, i, norbs, U, Up, Jh)
            hoppingTerms!(H, tarr, a, abit, i, norbs, N, state_list)
        end
    end

    return H
end


####################################################################################################
############ WITH SOC

function Hubbard_nfsz_multi_pbc(N::Int64,
    norbs::Int64,
    U::Float64,
    Up::Float64,
    Jh::Float64,
    tarr::Array{Tuple{Float64,Float64},1},
    state_list::Array{Int64,1}; lambda::Float64=0.0, choice::String="real")



    M = length(state_list)

    H = spzeros(Complex{Float64}, M, M)
    tarr = (N != 2) ? tarr : [0.5 .* i for i in tarr]

    for a in 1:M

        abit = digits(state_list[a], base=2, pad=2 * norbs * N)

        for i in 1:2*norbs:2*norbs*N-1

            #############  Interaction terms   ###############
            for m in 'a':'a'+norbs-1

                #U-term
                H[a, a] = H[a, a] + (U * (abit[ind(i, m).up] * abit[ind(i, m).down]))

                for n in 'a':'a'+norbs-1

                    if (m != n)

                        # (U'-J) and U'-terms
                        H[a, a] = H[a, a] + (0.5 * (Up - Jh) * ((abit[ind(i, m).up] * abit[ind(i, n).up]) + (abit[ind(i, m).down] * abit[ind(i, n).down])))
                        H[a, a] = H[a, a] + Up * abit[ind(i, m).up] * abit[ind(i, n).down]

                        # J-terms: spin-exchange and pair hopping
                        if (abit[ind(i, m).up] != abit[ind(i, n).up] && abit[ind(i, m).down] != abit[ind(i, n).down])

                            s_star = flip(flip(abit, ind(i, m).up, ind(i, n).up)[1], ind(i, m).down, ind(i, n).down)
                            b = findstate(s_star[2], M, state_list)
                            jw_str = (-1)^(abit[ind(i, m).up] + s_star[1][ind(i, n).up])

                            if (b > 0)
                                if (abit[ind(i, m).up] + abit[ind(i, m).down] == abit[ind(i, n).up] + abit[ind(i, n).down])
                                    H[a, b] = H[a, b] - 0.5 * jw_str * Jh #spin-exchange
                                else
                                    H[a, b] = H[a, b] + 0.5 * jw_str * Jh #pair-hop
                                end
                            end
                        end
                    end # end of "if(m!=n)"
                end  #for n in  'a':'a'+norbs-1
            end    #for m in 'a':'a'+norbs-1

            #############  Hopping terms   ###############

            ref_ind_site = i

            for t in tarr

                j = (i != 2 * norbs * N - (2 * norbs - 1)) ? (neigh(i, norbs)) : 1
                tm = t[1]
                tmn = t[2]

                for m in 'a':'a'+norbs-1

                    # Intra-orbital: i.e, tm-terms
                    # up-spin hopping

                    if (abit[ind(i, m).up] != abit[ind(j, m).up])
                        s_star = flip(abit, ind(i, m).up, ind(j, m).up)
                        b = findstate(s_star[2], M, state_list)
                        jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).up+1:ind(j, m).up-1) : (-1)^sum(abit[p] for p in ind(j, m).up+1:ind(i, m).up-1)

                        if (b > 0)
                            H[a, b] = H[a, b] - tm * jw_str
                        end
                    end

                    # down-spin hopping
                    if (abit[ind(i, m).down] != abit[ind(j, m).down])

                        s_star = flip(abit, ind(i, m).down, ind(j, m).down)
                        b = findstate(s_star[2], M, state_list)
                        jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).down+1:ind(j, m).down-1) : (-1)^sum(abit[p] for p in ind(j, m).down+1:ind(i, m).down-1)

                        if (b > 0)
                            H[a, b] = H[a, b] - tm * jw_str
                        end
                    end

                    # Inter-orbital: i.e, tmn-terms

                    for n in 'a':'a'+norbs-1

                        if (m != n)

                            # up-spin hopping
                            if (abit[ind(i, m).up] != abit[ind(j, n).up])
                                s_star = flip(abit, ind(i, m).up, ind(j, n).up)
                                b = findstate(s_star[2], M, state_list)
                                jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).up+1:ind(j, n).up-1) : (-1)^sum(abit[p] for p in ind(j, n).up+1:ind(i, m).up-1)

                                if (b > 0)
                                    H[a, b] = H[a, b] - tmn * jw_str
                                end
                            end

                            # down-spin hopping
                            if (abit[ind(i, m).down] != abit[ind(j, n).down])

                                s_star = flip(abit, ind(i, m).down, ind(j, n).down)
                                b = findstate(s_star[2], M, state_list)
                                jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).down+1:ind(j, n).down-1) : (-1)^sum(s_star[1][p] for p in ind(j, n).down+1:ind(i, m).down-1)

                                if (b > 0)
                                    H[a, b] = H[a, b] - tmn * jw_str
                                end
                            end
                        end #end of if(m!=n)
                    end  # for n in 'a':'a'+norbs-1 loop
                end #for m in 'a':'a'+norbs-1   loop

            end # #end of t loop

            if (lambda != 0.0)

                # SOC_x
                if (abit[ind(i, 'b').up] != abit[ind(i, 'c').down])

                    s_star = flip(abit, ind(i, 'b').up, ind(i, 'c').down)
                    jw_str = (-1)^sum(abit[p] for p in ind(i, 'b').up+1:ind(i, 'c').down-1)
                    b = findstate(s_star[2], M, state_list)

                    if (b > 0)
                        H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'c').down]) * lambda * jw_str * Lsign
                    end

                elseif (abit[ind(i, 'b').down] != abit[ind(i, 'c').up])

                    s_star = flip(abit, ind(i, 'b').down, ind(i, 'c').up)
                    jw_str = 1#(-1)^sum(abit[p] for p in ind(i,'b').down+1:ind(i,'c').up-1)
                    b = findstate(s_star[2], M, state_list)

                    if (b > 0)
                        H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'c').up]) * lambda * jw_str * Lsign
                    end

                end

                #SOC_y
                if (abit[ind(i, 'c').up] != abit[ind(i, 'a').down])

                    s_star = flip(abit, ind(i, 'c').up, ind(i, 'a').down)
                    jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').down+1:ind(i, 'c').up-1)
                    b = findstate(s_star[2], M, state_list)

                    if (b > 0)
                        H[a, b] = H[a, b] - 0.5 * lambda * jw_str * Lsign
                    end

                elseif (abit[ind(i, 'c').down] != abit[ind(i, 'a').up])

                    s_star = flip(abit, ind(i, 'c').down, ind(i, 'a').up)
                    jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'c').down-1)
                    b = findstate(s_star[2], M, state_list)

                    if (b > 0)
                        H[a, b] = H[a, b] + 0.5 * lambda * jw_str * Lsign
                    end
                end

                #SOC_z                   
                if (abit[ind(i, 'a').up] != abit[ind(i, 'b').up])

                    s_star = flip(abit, ind(i, 'a').up, ind(i, 'b').up)
                    jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'b').up-1)
                    b = findstate(s_star[2], M, state_list)

                    if (b > 0)
                        H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'b').up]) * lambda * jw_str * Lsign
                    end

                elseif (abit[ind(i, 'a').down] != abit[ind(i, 'b').down])

                    s_star = flip(abit, ind(i, 'a').down, ind(i, 'b').down)
                    jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').down+1:ind(i, 'b').down-1)
                    b = findstate(s_star[2], M, state_list)

                    if (b > 0)
                        H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'a').down]) * lambda * jw_str * Lsign
                    end
                end
            end  # end of if lambda !=0.0 loop                                                                                                          
        end    # for i in 1:2*norbs:2*norbs*N-1 loop
    end #for a in 1:M loop

    return H
end

#####################################################################################################################################################

function Hubbard_nfsz_multi_obc(N::Int64,
    norbs::Int64,
    U::Float64,
    Up::Float64,
    Jh::Float64,
    tarr::Array{Tuple{Float64,Float64},1},
    state_list::Array{Int64,1}; lambda::Float64=0.0, choice::String="real")



    M = length(state_list)

    H = spzeros(Complex{Float64}, M, M)
    tarr = (N != 2) ? tarr : [0.5 .* i for i in tarr]

    for a in 1:M

        abit = digits(state_list[a], base=2, pad=2 * norbs * N)

        for i in 1:2*norbs:2*norbs*N-1

            #############  Interaction terms   ###############
            for m in 'a':'a'+norbs-1

                #U-term
                H[a, a] = H[a, a] + (U * (abit[ind(i, m).up] * abit[ind(i, m).down]))

                for n in 'a':'a'+norbs-1

                    if (m != n)

                        # (U'-J) and U'-terms
                        H[a, a] = H[a, a] + (0.5 * (Up - Jh) * ((abit[ind(i, m).up] * abit[ind(i, n).up]) + (abit[ind(i, m).down] * abit[ind(i, n).down])))
                        H[a, a] = H[a, a] + Up * abit[ind(i, m).up] * abit[ind(i, n).down]

                        # J-terms: spin-exchange and pair hopping
                        if (abit[ind(i, m).up] != abit[ind(i, n).up] && abit[ind(i, m).down] != abit[ind(i, n).down])

                            s_star = flip(flip(abit, ind(i, m).up, ind(i, n).up)[1], ind(i, m).down, ind(i, n).down)
                            b = findstate(s_star[2], M, state_list)
                            jw_str = (-1)^(abit[ind(i, m).up] + s_star[1][ind(i, n).up])

                            if (b > 0)
                                if (abit[ind(i, m).up] + abit[ind(i, m).down] == abit[ind(i, n).up] + abit[ind(i, n).down])
                                    H[a, b] = H[a, b] - 0.5 * jw_str * Jh #spin-exchange
                                else
                                    H[a, b] = H[a, b] + 0.5 * jw_str * Jh #pair-hop
                                end
                            end
                        end
                    end # end of "if(m!=n)"
                end  #for n in  'a':'a'+norbs-1
            end    #for m in 'a':'a'+norbs-1

            #############  Hopping terms   ###############

            ref_ind_site = i

            for t in tarr

                j = (i != 2 * norbs * N - (2 * norbs - 1)) ? (neigh(i, norbs)) : @goto here

                tm = t[1]
                tmn = t[2]

                for m in 'a':'a'+norbs-1

                    # Intra-orbital: i.e, tm-terms
                    # up-spin hopping

                    if (abit[ind(i, m).up] != abit[ind(j, m).up])
                        s_star = flip(abit, ind(i, m).up, ind(j, m).up)
                        b = findstate(s_star[2], M, state_list)
                        jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).up+1:ind(j, m).up-1) : (-1)^sum(abit[p] for p in ind(j, m).up+1:ind(i, m).up-1)

                        if (b > 0)
                            H[a, b] = H[a, b] - tm * jw_str
                        end
                    end

                    # down-spin hopping
                    if (abit[ind(i, m).down] != abit[ind(j, m).down])

                        s_star = flip(abit, ind(i, m).down, ind(j, m).down)
                        b = findstate(s_star[2], M, state_list)
                        jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).down+1:ind(j, m).down-1) : (-1)^sum(abit[p] for p in ind(j, m).down+1:ind(i, m).down-1)

                        if (b > 0)
                            H[a, b] = H[a, b] - tm * jw_str
                        end
                    end

                    # Inter-orbital: i.e, tmn-terms

                    for n in 'a':'a'+norbs-1

                        if (m != n)

                            # up-spin hopping
                            if (abit[ind(i, m).up] != abit[ind(j, n).up])
                                s_star = flip(abit, ind(i, m).up, ind(j, n).up)
                                b = findstate(s_star[2], M, state_list)
                                jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).up+1:ind(j, n).up-1) : (-1)^sum(abit[p] for p in ind(j, n).up+1:ind(i, m).up-1)

                                if (b > 0)
                                    H[a, b] = H[a, b] - tmn * jw_str
                                end
                            end

                            # down-spin hopping
                            if (abit[ind(i, m).down] != abit[ind(j, n).down])

                                s_star = flip(abit, ind(i, m).down, ind(j, n).down)
                                b = findstate(s_star[2], M, state_list)
                                jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).down+1:ind(j, n).down-1) : (-1)^sum(s_star[1][p] for p in ind(j, n).down+1:ind(i, m).down-1)

                                if (b > 0)
                                    H[a, b] = H[a, b] - tmn * jw_str
                                end
                            end
                        end #end of if(m!=n)
                    end  # for n in 'a':'a'+norbs-1 loop
                end #for m in 'a':'a'+norbs-1   loop

            end # #end of t loop
            @label here

            if (lambda != 0.0)

                # SOC_x
                if (abit[ind(i, 'b').up] != abit[ind(i, 'c').down])

                    s_star = flip(abit, ind(i, 'b').up, ind(i, 'c').down)
                    jw_str = (-1)^sum(abit[p] for p in ind(i, 'b').up+1:ind(i, 'c').down-1)
                    b = findstate(s_star[2], M, state_list)

                    if (b > 0)
                        H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'c').down]) * lambda * jw_str * Lsign
                    end

                elseif (abit[ind(i, 'b').down] != abit[ind(i, 'c').up])

                    s_star = flip(abit, ind(i, 'b').down, ind(i, 'c').up)
                    jw_str = 1#(-1)^sum(abit[p] for p in ind(i,'b').down+1:ind(i,'c').up-1)
                    b = findstate(s_star[2], M, state_list)

                    if (b > 0)
                        H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'c').up]) * lambda * jw_str * Lsign
                    end

                end

                #SOC_y
                if (abit[ind(i, 'c').up] != abit[ind(i, 'a').down])

                    s_star = flip(abit, ind(i, 'c').up, ind(i, 'a').down)
                    jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').down+1:ind(i, 'c').up-1)
                    b = findstate(s_star[2], M, state_list)

                    if (b > 0)
                        H[a, b] = H[a, b] - 0.5 * lambda * jw_str * Lsign
                    end

                elseif (abit[ind(i, 'c').down] != abit[ind(i, 'a').up])

                    s_star = flip(abit, ind(i, 'c').down, ind(i, 'a').up)
                    jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'c').down-1)
                    b = findstate(s_star[2], M, state_list)

                    if (b > 0)
                        H[a, b] = H[a, b] + 0.5 * lambda * jw_str * Lsign
                    end
                end

                #SOC_z                   
                if (abit[ind(i, 'a').up] != abit[ind(i, 'b').up])

                    s_star = flip(abit, ind(i, 'a').up, ind(i, 'b').up)
                    jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'b').up-1)
                    b = findstate(s_star[2], M, state_list)

                    if (b > 0)
                        H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'b').up]) * lambda * jw_str * Lsign
                    end

                elseif (abit[ind(i, 'a').down] != abit[ind(i, 'b').down])

                    s_star = flip(abit, ind(i, 'a').down, ind(i, 'b').down)
                    jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').down+1:ind(i, 'b').down-1)
                    b = findstate(s_star[2], M, state_list)

                    if (b > 0)
                        H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'a').down]) * lambda * jw_str * Lsign
                    end
                end
            end  # end of if lambda !=0.0 loop                                                                                                          
        end    # for i in 1:2*norbs:2*norbs*N-1 loop
    end #for a in 1:M loop

    return H
end


#####################################################################################################################################################

# (2b) Hamiltonian for momentum conservation (with additional symmetries like nf/nfsz)

function Hubbard_tr_multi(N::Int64,
    norbs::Int64,
    U::Float64,
    Up::Float64,
    Jh::Float64,
    tarr::Array{Tuple{Float64,Float64},1},
    k::Int64,
    nf::Int64,
    state_listinfo::Tuple{Array{Int64,1},Array{Float64,1}})



    rep_states = state_listinfo[1]
    Nconst = state_listinfo[2]

    M = length(rep_states)
    H = Complex.(spzeros(M, M))
    trfac = exp(-im * 2 * pi * k / N)
    tarr = (N != 2) ? tarr : [0.5 .* i for i in tarr]

    # For even fermion sectors, an additional sign is carried depending upon the respective state and its relation with its representative
    if (nf % 2 == 0)

        for a in 1:M

            abit = digits(rep_states[a], base=2, pad=2 * norbs * N)

            for i in 1:2*norbs:2*norbs*N-1

                #############  Interaction terms   ###############
                for m in 'a':'a'+norbs-1

                    #U-term
                    H[a, a] = H[a, a] + (U * (abit[ind(i, m).up] * abit[ind(i, m).down]))

                    for n in 'a':'a'+norbs-1

                        if (m != n)

                            # (U'-J) and U'-terms
                            H[a, a] = H[a, a] + (0.5 * (Up - Jh) * ((abit[ind(i, m).up] * abit[ind(i, n).up]) + (abit[ind(i, m).down] * abit[ind(i, n).down])))
                            H[a, a] = H[a, a] + Up * abit[ind(i, m).up] * abit[ind(i, n).down]

                            # J-terms: spin-exchange and pair hopping
                            if (abit[ind(i, m).up] != abit[ind(i, n).up] && abit[ind(i, m).down] != abit[ind(i, n).down])

                                s_star = flip(flip(abit, ind(i, m).up, ind(i, n).up)[1], ind(i, m).down, ind(i, n).down)
                                b_info = representative_even(s_star[2], s_star[1], N, norbs)
                                b = findstate(b_info[1], M, rep_states)
                                jw_str = (-1)^(abit[ind(i, m).up] + s_star[1][ind(i, n).up])

                                if (b > 0)
                                    if (abit[ind(i, m).up] + abit[ind(i, m).down] == abit[ind(i, n).up] + abit[ind(i, n).down])
                                        H[a, b] = H[a, b] - 0.5 * Jh * jw_str * b_info[3] * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    else
                                        H[a, b] = H[a, b] + 0.5 * Jh * jw_str * b_info[3] * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    end
                                end
                            end
                        end # end of "if(m!=n)"
                    end  #for n in  'a':'a'+norbs-1
                end    #for m in 'a':'a'+norbs-1

                #############  Hopping terms   ###############



                for t in tarr


                    j = (i != 2 * norbs * N - (2 * norbs - 1)) ? (neigh(i, norbs)) : 1

                    tm = t[1]
                    tmm = t[2]

                    for m in 'a':'a'+norbs-1

                        # Intra-orbital: i.e, tm-terms
                        # up-spin hopping
                        if (abit[ind(i, m).up] != abit[ind(j, m).up])
                            s_star = flip(abit, ind(i, m).up, ind(j, m).up)
                            b_info = representative_even(s_star[2], s_star[1], N, norbs)
                            b = findstate(b_info[1], M, rep_states)
                            jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).up+1:ind(j, m).up-1) : (-1)^sum(abit[p] for p in ind(j, m).up+1:ind(i, m).up-1)

                            if (b > 0)
                                H[a, b] = H[a, b] - tm * jw_str * b_info[3] * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                            end
                        end

                        # down-spin hopping
                        if (abit[ind(i, m).down] != abit[ind(j, m).down])

                            s_star = flip(abit, ind(i, m).down, ind(j, m).down)
                            b_info = representative_even(s_star[2], s_star[1], N, norbs)
                            b = findstate(b_info[1], M, rep_states)
                            jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).down+1:ind(j, m).down-1) : (-1)^sum(abit[p] for p in ind(j, m).down+1:ind(i, m).down-1)

                            if (b > 0)
                                H[a, b] = H[a, b] - tm * jw_str * b_info[3] * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                            end
                        end

                        # Inter-orbital: i.e, tmn-terms

                        for n in 'a':'a'+norbs-1

                            if (m != n)

                                # up-spin hopping
                                if (abit[ind(i, m).up] != abit[ind(j, n).up])
                                    s_star = flip(abit, ind(i, m).up, ind(j, n).up)
                                    b_info = representative_even(s_star[2], s_star[1], N, norbs)
                                    b = findstate(b_info[1], M, rep_states)
                                    jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).up+1:ind(j, n).up-1) : (-1)^sum(abit[p] for p in ind(j, n).up+1:ind(i, m).up-1)

                                    if (b > 0)
                                        H[a, b] = H[a, b] - tmm * jw_str * b_info[3] * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    end
                                end

                                # down-spin hopping
                                if (abit[ind(i, m).down] != abit[ind(j, n).down])

                                    s_star = flip(abit, ind(i, m).down, ind(j, n).down)
                                    b_info = representative_even(s_star[2], s_star[1], N, norbs)
                                    b = findstate(b_info[1], M, rep_states)
                                    jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).down+1:ind(j, n).down-1) : (-1)^sum(abit[p] for p in ind(j, n).down+1:ind(i, m).down-1)

                                    if (b > 0)
                                        H[a, b] = H[a, b] - tmm * jw_str * b_info[3] * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    end
                                end
                            end #end of if(m!=n)
                        end  # for n in 'a':'a'+norbs-1 loop
                    end #for m in 'a':'a'+norbs-1   loop

                end # #end of t loop

            end    # for i in 1:2*norbs:2*norbs*N-1 loop
        end #for a in 1:M loop

    else

        for a in 1:M

            abit = digits(rep_states[a], base=2, pad=2 * norbs * N)

            for i in 1:2*norbs:2*norbs*N-1

                for m in 'a':'a'+norbs-1

                    #U-term
                    H[a, a] = H[a, a] + (U * (abit[ind(i, m).up] * abit[ind(i, m).down]))

                    for n in 'a':'a'+norbs-1

                        if (m != n)

                            H[a, a] = H[a, a] + (0.5 * (Up - Jh) * ((abit[ind(i, m).up] * abit[ind(i, n).up]) + (abit[ind(i, m).down] * abit[ind(i, n).down])))
                            H[a, a] = H[a, a] + Up * abit[ind(i, m).up] * abit[ind(i, n).down]

                            if (abit[ind(i, m).up] != abit[ind(i, n).up] && abit[ind(i, m).down] != abit[ind(i, n).down])

                                s_star = flip(flip(abit, ind(i, m).up, ind(i, n).up)[1], ind(i, m).down, ind(i, n).down)
                                b_info = representative_odd(s_star[2], N, norbs)
                                b = findstate(b_info[1], M, rep_states)
                                jw_str = (-1)^(abit[ind(i, m).up] + s_star[1][ind(i, n).up])
                                if (b > 0)
                                    if (abit[ind(i, m).up] + abit[ind(i, m).down] == abit[ind(i, n).up] + abit[ind(i, n).down])
                                        H[a, b] = H[a, b] - 0.5 * Jh * jw_str * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    else
                                        H[a, b] = H[a, b] + 0.5 * Jh * jw_str * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    end
                                end
                            end
                        end # end of "if(m!=n)"
                    end  #for n in  'a':'a'+norbs-1
                end    #for m in 'a':'a'+norbs-1

                #t-terms

                for t in tarr

                    j = (i != 2 * norbs * N - (2 * norbs - 1)) ? (neigh(i, norbs)) : 1

                    tm = t[1]
                    tmm = t[2]

                    for m in 'a':'a'+norbs-1

                        if (abit[ind(i, m).up] != abit[ind(j, m).up]) #up-spin
                            s_star = flip(abit, ind(i, m).up, ind(j, m).up)
                            b_info = representative_odd(s_star[2], N, norbs)
                            b = findstate(b_info[1], M, rep_states)
                            jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).up+1:ind(j, m).up-1) : (-1)^sum(abit[p] for p in ind(j, m).up+1:ind(i, m).up-1)

                            if (b > 0)
                                H[a, b] = H[a, b] - tm * jw_str * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                            end
                        end

                        if (abit[ind(i, m).down] != abit[ind(j, m).down])  #down-spin

                            s_star = flip(abit, ind(i, m).down, ind(j, m).down)
                            b_info = representative_odd(s_star[2], N, norbs)
                            b = findstate(b_info[1], M, rep_states)
                            jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).down+1:ind(j, m).down-1) : (-1)^sum(abit[p] for p in ind(j, m).down+1:ind(i, m).down-1)

                            if (b > 0)
                                H[a, b] = H[a, b] - tm * jw_str * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                            end
                        end


                        for n in 'a':'a'+norbs-1

                            if (m != n)

                                if (abit[ind(i, m).up] != abit[ind(j, n).up]) #up-spin
                                    s_star = flip(abit, ind(i, m).up, ind(j, n).up)
                                    b_info = representative_odd(s_star[2], N, norbs)
                                    b = findstate(b_info[1], M, rep_states)
                                    jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).up+1:ind(j, n).up-1) : (-1)^sum(abit[p] for p in ind(j, n).up+1:ind(i, m).up-1)

                                    if (b > 0)
                                        H[a, b] = H[a, b] - tmm * jw_str * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    end
                                end

                                if (abit[ind(i, m).down] != abit[ind(j, n).down]) #down-spin

                                    s_star = flip(abit, ind(i, m).down, ind(j, n).down)
                                    b_info = representative_odd(s_star[2], N, norbs)
                                    b = findstate(b_info[1], M, rep_states)
                                    jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).down+1:ind(j, n).down-1) : (-1)^sum(abit[p] for p in ind(j, n).down+1:ind(i, m).down-1)

                                    if (b > 0)
                                        H[a, b] = H[a, b] - tmm * jw_str * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    end
                                end
                            end #end of if(m!=n)
                        end  # for n in 'a':'a'+norbs-1 loop
                    end #for m in 'a':'a'+norbs-1   loop

                end # #end of t loop

            end    # for i in 1:2*norbs:2*norbs*N-1 loop
        end #for a in 1:M loop
    end

    return H
end

##############################################################################################################
######## WITH SOC

function Hubbard_tr_multi(N::Int64,
    norbs::Int64,
    U::Float64,
    Up::Float64,
    Jh::Float64,
    tarr::Array{Tuple{Float64,Float64},1},
    k::Int64,
    nf::Int64,
    state_listinfo::Tuple{Array{Int64,1},Array{Float64,1}}; lambda::Float64=0.0, choice::String="real")



    rep_states = state_listinfo[1]
    Nconst = state_listinfo[2]

    M = length(rep_states)
    H = Complex.(spzeros(M, M))
    trfac = exp(-im * 2 * pi * k / N)
    tarr = (N != 2) ? tarr : [0.5 .* i for i in tarr]

    # For even fermion sectors, an additional sign is carried depending upon the respective state and its relation with its representative
    if (nf % 2 == 0)

        for a in 1:M

            abit = digits(rep_states[a], base=2, pad=2 * norbs * N)

            for i in 1:2*norbs:2*norbs*N-1

                #############  Interaction terms   ###############
                for m in 'a':'a'+norbs-1

                    #U-term
                    H[a, a] = H[a, a] + (U * (abit[ind(i, m).up] * abit[ind(i, m).down]))

                    for n in 'a':'a'+norbs-1

                        if (m != n)

                            # (U'-J) and U'-terms
                            H[a, a] = H[a, a] + (0.5 * (Up - Jh) * ((abit[ind(i, m).up] * abit[ind(i, n).up]) + (abit[ind(i, m).down] * abit[ind(i, n).down])))
                            H[a, a] = H[a, a] + Up * abit[ind(i, m).up] * abit[ind(i, n).down]

                            # J-terms: spin-exchange and pair hopping
                            if (abit[ind(i, m).up] != abit[ind(i, n).up] && abit[ind(i, m).down] != abit[ind(i, n).down])

                                s_star = flip(flip(abit, ind(i, m).up, ind(i, n).up)[1], ind(i, m).down, ind(i, n).down)
                                b_info = representative_even(s_star[2], s_star[1], N, norbs)
                                b = findstate(b_info[1], M, rep_states)
                                jw_str = (-1)^(abit[ind(i, m).up] + s_star[1][ind(i, n).up])

                                if (b > 0)
                                    if (abit[ind(i, m).up] + abit[ind(i, m).down] == abit[ind(i, n).up] + abit[ind(i, n).down])
                                        H[a, b] = H[a, b] - 0.5 * Jh * jw_str * b_info[3] * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    else
                                        H[a, b] = H[a, b] + 0.5 * Jh * jw_str * b_info[3] * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    end
                                end
                            end
                        end # end of "if(m!=n)"
                    end  #for n in  'a':'a'+norbs-1
                end    #for m in 'a':'a'+norbs-1

                #############  Hopping terms   ###############



                for t in tarr


                    j = (i != 2 * norbs * N - (2 * norbs - 1)) ? (neigh(i, norbs)) : 1

                    tm = t[1]
                    tmm = t[2]

                    for m in 'a':'a'+norbs-1

                        # Intra-orbital: i.e, tm-terms
                        # up-spin hopping
                        if (abit[ind(i, m).up] != abit[ind(j, m).up])
                            s_star = flip(abit, ind(i, m).up, ind(j, m).up)
                            b_info = representative_even(s_star[2], s_star[1], N, norbs)
                            b = findstate(b_info[1], M, rep_states)
                            jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).up+1:ind(j, m).up-1) : (-1)^sum(abit[p] for p in ind(j, m).up+1:ind(i, m).up-1)

                            if (b > 0)
                                H[a, b] = H[a, b] - tm * jw_str * b_info[3] * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                            end
                        end

                        # down-spin hopping
                        if (abit[ind(i, m).down] != abit[ind(j, m).down])

                            s_star = flip(abit, ind(i, m).down, ind(j, m).down)
                            b_info = representative_even(s_star[2], s_star[1], N, norbs)
                            b = findstate(b_info[1], M, rep_states)
                            jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).down+1:ind(j, m).down-1) : (-1)^sum(abit[p] for p in ind(j, m).down+1:ind(i, m).down-1)

                            if (b > 0)
                                H[a, b] = H[a, b] - tm * jw_str * b_info[3] * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                            end
                        end

                        # Inter-orbital: i.e, tmn-terms

                        for n in 'a':'a'+norbs-1

                            if (m != n)

                                # up-spin hopping
                                if (abit[ind(i, m).up] != abit[ind(j, n).up])
                                    s_star = flip(abit, ind(i, m).up, ind(j, n).up)
                                    b_info = representative_even(s_star[2], s_star[1], N, norbs)
                                    b = findstate(b_info[1], M, rep_states)
                                    jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).up+1:ind(j, n).up-1) : (-1)^sum(abit[p] for p in ind(j, n).up+1:ind(i, m).up-1)

                                    if (b > 0)
                                        H[a, b] = H[a, b] - tmm * jw_str * b_info[3] * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    end
                                end

                                # down-spin hopping
                                if (abit[ind(i, m).down] != abit[ind(j, n).down])

                                    s_star = flip(abit, ind(i, m).down, ind(j, n).down)
                                    b_info = representative_even(s_star[2], s_star[1], N, norbs)
                                    b = findstate(b_info[1], M, rep_states)
                                    jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).down+1:ind(j, n).down-1) : (-1)^sum(abit[p] for p in ind(j, n).down+1:ind(i, m).down-1)

                                    if (b > 0)
                                        H[a, b] = H[a, b] - tmm * jw_str * b_info[3] * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    end
                                end
                            end #end of if(m!=n)
                        end  # for n in 'a':'a'+norbs-1 loop
                    end #for m in 'a':'a'+norbs-1   loop

                end # #end of t loop

                if (lambda != 0.0)

                    # SOC_x
                    if (abit[ind(i, 'b').up] != abit[ind(i, 'c').down])

                        s_star = flip(abit, ind(i, 'b').up, ind(i, 'c').down)
                        jw_str = (-1)^sum(abit[p] for p in ind(i, 'b').up+1:ind(i, 'c').down-1)
                        b_info = representative_even(s_star[2], N, norbs)
                        b = findstate(b_info[1], M, state_list)

                        if (b > 0)
                            H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'c').down]) * lambda * jw_str * Lsign * b_info[3] * (Nconst[b] / Nconst[a]) * exp((-im * 2 * pi * k * b_info[2]) / N)
                        end

                    elseif (abit[ind(i, 'b').down] != abit[ind(i, 'c').up])

                        s_star = flip(abit, ind(i, 'b').down, ind(i, 'c').up)
                        jw_str = 1#(-1)^sum(abit[p] for p in ind(i,'b').down+1:ind(i,'c').up-1)
                        b_info = representative_even(s_star[2], N, norbs)
                        b = findstate(b_info[1], M, state_list)

                        if (b > 0)
                            H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'c').up]) * lambda * jw_str * Lsign * b_info[3] * (Nconst[b] / Nconst[a]) * exp((-im * 2 * pi * k * b_info[2]) / N)
                        end

                    end

                    #SOC_y
                    if (abit[ind(i, 'c').up] != abit[ind(i, 'a').down])

                        s_star = flip(abit, ind(i, 'c').up, ind(i, 'a').down)
                        jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').down+1:ind(i, 'c').up-1)
                        b_info = representative_even(s_star[2], N, norbs)
                        b = findstate(b_info[1], M, state_list)

                        if (b > 0)
                            H[a, b] = H[a, b] - 0.5 * lambda * jw_str * Lsign * b_info[3] * (Nconst[b] / Nconst[a]) * exp((-im * 2 * pi * k * b_info[2]) / N)
                        end

                    elseif (abit[ind(i, 'c').down] != abit[ind(i, 'a').up])

                        s_star = flip(abit, ind(i, 'c').down, ind(i, 'a').up)
                        jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'c').down-1)
                        b_info = representative_even(s_star[2], N, norbs)
                        b = findstate(b_info[1], M, state_list)

                        if (b > 0)
                            H[a, b] = H[a, b] + 0.5 * lambda * jw_str * Lsign * b_info[3] * (Nconst[b] / Nconst[a]) * exp((-im * 2 * pi * k * b_info[2]) / N)
                        end
                    end

                    #SOC_z                   
                    if (abit[ind(i, 'a').up] != abit[ind(i, 'b').up])

                        s_star = flip(abit, ind(i, 'a').up, ind(i, 'b').up)
                        jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'b').up-1)
                        b_info = representative_even(s_star[2], N, norbs)
                        b = findstate(b_info[1], M, state_list)

                        if (b > 0)
                            H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'b').up]) * lambda * jw_str * Lsign * b_info[3] * (Nconst[b] / Nconst[a]) * exp((-im * 2 * pi * k * b_info[2]) / N)
                        end

                    elseif (abit[ind(i, 'a').down] != abit[ind(i, 'b').down])

                        s_star = flip(abit, ind(i, 'a').down, ind(i, 'b').down)
                        jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').down+1:ind(i, 'b').down-1)
                        b_info = representative_even(s_star[2], N, norbs)
                        b = findstate(b_info[1], M, state_list)

                        if (b > 0)
                            H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'a').down]) * lambda * jw_str * Lsign * b_info[3] * (Nconst[b] / Nconst[a]) * exp((-im * 2 * pi * k * b_info[2]) / N)
                        end
                    end
                end  # end of if lambda !=0.0 loop
            end    # for i in 1:2*norbs:2*norbs*N-1 loop
        end #for a in 1:M loop

    else

        for a in 1:M

            abit = digits(rep_states[a], base=2, pad=2 * norbs * N)

            for i in 1:2*norbs:2*norbs*N-1

                for m in 'a':'a'+norbs-1

                    #U-term
                    H[a, a] = H[a, a] + (U * (abit[ind(i, m).up] * abit[ind(i, m).down]))

                    for n in 'a':'a'+norbs-1

                        if (m != n)

                            H[a, a] = H[a, a] + (0.5 * (Up - Jh) * ((abit[ind(i, m).up] * abit[ind(i, n).up]) + (abit[ind(i, m).down] * abit[ind(i, n).down])))
                            H[a, a] = H[a, a] + Up * abit[ind(i, m).up] * abit[ind(i, n).down]

                            if (abit[ind(i, m).up] != abit[ind(i, n).up] && abit[ind(i, m).down] != abit[ind(i, n).down])

                                s_star = flip(flip(abit, ind(i, m).up, ind(i, n).up)[1], ind(i, m).down, ind(i, n).down)
                                b_info = representative_odd(s_star[2], N, norbs)
                                b = findstate(b_info[1], M, rep_states)
                                jw_str = (-1)^(abit[ind(i, m).up] + s_star[1][ind(i, n).up])
                                if (b > 0)
                                    if (abit[ind(i, m).up] + abit[ind(i, m).down] == abit[ind(i, n).up] + abit[ind(i, n).down])
                                        H[a, b] = H[a, b] - 0.5 * Jh * jw_str * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    else
                                        H[a, b] = H[a, b] + 0.5 * Jh * jw_str * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    end
                                end
                            end
                        end # end of "if(m!=n)"
                    end  #for n in  'a':'a'+norbs-1
                end    #for m in 'a':'a'+norbs-1

                #t-terms

                for t in tarr

                    j = (i != 2 * norbs * N - (2 * norbs - 1)) ? (neigh(i, norbs)) : 1

                    tm = t[1]
                    tmm = t[2]

                    for m in 'a':'a'+norbs-1

                        if (abit[ind(i, m).up] != abit[ind(j, m).up]) #up-spin
                            s_star = flip(abit, ind(i, m).up, ind(j, m).up)
                            b_info = representative_odd(s_star[2], N, norbs)
                            b = findstate(b_info[1], M, rep_states)
                            jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).up+1:ind(j, m).up-1) : (-1)^sum(abit[p] for p in ind(j, m).up+1:ind(i, m).up-1)

                            if (b > 0)
                                H[a, b] = H[a, b] - tm * jw_str * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                            end
                        end

                        if (abit[ind(i, m).down] != abit[ind(j, m).down])  #down-spin

                            s_star = flip(abit, ind(i, m).down, ind(j, m).down)
                            b_info = representative_odd(s_star[2], N, norbs)
                            b = findstate(b_info[1], M, rep_states)
                            jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).down+1:ind(j, m).down-1) : (-1)^sum(abit[p] for p in ind(j, m).down+1:ind(i, m).down-1)

                            if (b > 0)
                                H[a, b] = H[a, b] - tm * jw_str * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                            end
                        end


                        for n in 'a':'a'+norbs-1

                            if (m != n)

                                if (abit[ind(i, m).up] != abit[ind(j, n).up]) #up-spin
                                    s_star = flip(abit, ind(i, m).up, ind(j, n).up)
                                    b_info = representative_odd(s_star[2], N, norbs)
                                    b = findstate(b_info[1], M, rep_states)
                                    jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).up+1:ind(j, n).up-1) : (-1)^sum(abit[p] for p in ind(j, n).up+1:ind(i, m).up-1)

                                    if (b > 0)
                                        H[a, b] = H[a, b] - tmm * jw_str * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    end
                                end

                                if (abit[ind(i, m).down] != abit[ind(j, n).down]) #down-spin

                                    s_star = flip(abit, ind(i, m).down, ind(j, n).down)
                                    b_info = representative_odd(s_star[2], N, norbs)
                                    b = findstate(b_info[1], M, rep_states)
                                    jw_str = i < j ? (-1)^sum(abit[p] for p in ind(i, m).down+1:ind(j, n).down-1) : (-1)^sum(abit[p] for p in ind(j, n).down+1:ind(i, m).down-1)

                                    if (b > 0)
                                        H[a, b] = H[a, b] - tmm * jw_str * (Nconst[b] / Nconst[a]) * (trfac^b_info[2])
                                    end
                                end
                            end #end of if(m!=n)
                        end  # for n in 'a':'a'+norbs-1 loop
                    end #for m in 'a':'a'+norbs-1   loop

                end # #end of t loop

                if (lambda != 0.0)

                    # SOC_x
                    if (abit[ind(i, 'b').up] != abit[ind(i, 'c').down])

                        s_star = flip(abit, ind(i, 'b').up, ind(i, 'c').down)
                        jw_str = (-1)^sum(abit[p] for p in ind(i, 'b').up+1:ind(i, 'c').down-1)
                        b_info = representative_odd(s_star[2], N, norbs)
                        b = findstate(b_info[1], M, rep_states)

                        if (b > 0)
                            H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'c').down]) * lambda * jw_str * Lsign * (Nconst[b] / Nconst[a]) * exp((-im * 2 * pi * k * b_info[2]) / N)
                        end

                    elseif (abit[ind(i, 'b').down] != abit[ind(i, 'c').up])

                        s_star = flip(abit, ind(i, 'b').down, ind(i, 'c').up)
                        jw_str = 1#(-1)^sum(abit[p] for p in ind(i,'b').down+1:ind(i,'c').up-1)
                        b_info = representative_odd(s_star[2], N, norbs)
                        b = findstate(b_info[1], M, rep_states)

                        if (b > 0)
                            H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'c').up]) * lambda * jw_str * Lsign * (Nconst[b] / Nconst[a]) * exp((-im * 2 * pi * k * b_info[2]) / N)
                        end

                    end

                    #SOC_y
                    if (abit[ind(i, 'c').up] != abit[ind(i, 'a').down])

                        s_star = flip(abit, ind(i, 'c').up, ind(i, 'a').down)
                        jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').down+1:ind(i, 'c').up-1)
                        b_info = representative_odd(s_star[2], N, norbs)
                        b = findstate(b_info[1], M, rep_states)

                        if (b > 0)
                            H[a, b] = H[a, b] - 0.5 * lambda * jw_str * Lsign * (Nconst[b] / Nconst[a]) * exp((-im * 2 * pi * k * b_info[2]) / N)
                        end

                    elseif (abit[ind(i, 'c').down] != abit[ind(i, 'a').up])

                        s_star = flip(abit, ind(i, 'c').down, ind(i, 'a').up)
                        jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'c').down-1)
                        b_info = representative_odd(s_star[2], N, norbs)
                        b = findstate(b_info[1], M, rep_states)

                        if (b > 0)
                            H[a, b] = H[a, b] + 0.5 * lambda * jw_str * Lsign * (Nconst[b] / Nconst[a]) * exp((-im * 2 * pi * k * b_info[2]) / N)
                        end
                    end

                    #SOC_z                   
                    if (abit[ind(i, 'a').up] != abit[ind(i, 'b').up])

                        s_star = flip(abit, ind(i, 'a').up, ind(i, 'b').up)
                        jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').up+1:ind(i, 'b').up-1)
                        b_info = representative_odd(s_star[2], N, norbs)
                        b = findstate(b_info[1], M, rep_states)

                        if (b > 0)
                            H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'b').up]) * lambda * jw_str * Lsign * (Nconst[b] / Nconst[a]) * exp((-im * 2 * pi * k * b_info[2]) / N)
                        end

                    elseif (abit[ind(i, 'a').down] != abit[ind(i, 'b').down])

                        s_star = flip(abit, ind(i, 'a').down, ind(i, 'b').down)
                        jw_str = (-1)^sum(abit[p] for p in ind(i, 'a').down+1:ind(i, 'b').down-1)
                        b_info = representative_odd(s_star[2], N, norbs)
                        b = findstate(b_info[1], M, rep_states)

                        if (b > 0)
                            H[a, b] = H[a, b] + 0.5 * im * ((-1)^abit[ind(i, 'a').down]) * lambda * jw_str * Lsign * (Nconst[b] / Nconst[a]) * exp((-im * 2 * pi * k * b_info[2]) / N)
                        end
                    end
                end  # end of if lambda !=0.0 loop
            end    # for i in 1:2*norbs:2*norbs*N-1 loop
        end #for a in 1:M loop
    end

    return H
end
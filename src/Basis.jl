"""
Notes: Possibly change 
Output of ind to return a struct
"""
dict = Dict(j => i for (j, i) in zip('a':'e', 1:5))
ind(i, orbname) = (up=i + (2 * (dict[orbname] - 1)), down=i + (2 * (dict[orbname] - 1)) + 1)
neigh(i, norbs) = i + (2 * norbs)

"""
--------------------------------------------------------------------------------
Basis manipulation functions                                 
--------------------------------------------------------------------------------
"""

"""
    function flip(a::Vector{Int64}, i::Int64, j::Int64)
XOR the i,j index of given array with 1,
without repeated XOR in case of i=j.
"""
function flip(a::Vector{Int64}, i::Int64, j::Int64)

    a_flip = copy(a)
    if i == j
        a_flip[i] = xor(a_flip[i], 1)
    else
        a_flip[i] = xor(a_flip[i], 1)
        a_flip[j] = xor(a_flip[j], 1)
    end
    b = zero(Int64)
    for i in eachindex(a_flip)
        b += a_flip[i] * 2^(i - 1)
    end
    return a_flip, b
end

"""
    function flip(a::Int64, indexes::Vararg{Int64,N}) where {N}
Flip the states at each of the indexes.
In the case of overlap only flip the state once.
"""
function flip(a::Int64, indexes::Vararg{Int64,N}) where {N}
    mask = mapreduce((x) -> 1 << (x - 1), |, indexes, init=0)
    return xor(a, mask)
end

"""
    function cyclebits(t::Int64, N::Int64, norbs::Int64)
Perform translation (used in momentum-conservation)
"""
function cyclebits(t::Int64, N::Int64, norbs::Int64)
    padIndex = 0
    i = t
    while !(iszero(i))
        i >>= 1
        padIndex += 1
    end
    pad = max(2 * norbs * N, padIndex)
    shift = 2 * norbs
    return t >> shift | (t << (pad - shift) & (1 << pad - 1))
end

"""
    function findstate(s_star::Int64, state_list::Vector{Int64})
Finds the location of a state in given state list
It is a function common to Hamiltonians with and without tr-invariance
"""
function findstate(s_star::Int64, state_list::Vector{Int64})
    result = searchsorted(state_list, s_star)
    if isempty(result)
        return -1
    else
        return first(result)
    end
end

"""
--------------------------------------------------------------------------------
checkstate functions (used in momentum-conservation):
Check if a state can be a representative and return its periodicity
--------------------------------------------------------------------------------
"""

function checkstate(s::Int64, N::Int64, norbs::Int64, FaTerm::F) where {F<:Function}
    t = copy(s)
    states = Int64[]
    Fa = zero(ComplexF64)

    for i in 1:N
        t = cyclebits(t, N, norbs)
        if (t < s)
            @goto exitEarly
        elseif (t == s)
            Fa += FaTerm(i)
        elseif !(t in states)
            push!(states, t)
        end
    end

    if (abs2(Fa) > 10^-10)
        Da = length(states)
        return sqrt(abs2(Fa) * (Da + 1))
    else
        @label exitEarly
        return -one(Float64)
    end
end

function checkstate_even(s::Int64, N::Int64, k::Int64, norbs::Int64)
    function FaTerm(s, N, k, norbs)
        sbit = digits(s, base=2, pad=2 * norbs * N)
        return (i) -> (-1)^sum(sbit[1:2*norbs*i]) * exp(-im * ((2 * pi * k * i) / N))
    end
    return checkstate(s, N, norbs, FaTerm(s, N, k, norbs))
end

function checkstate_odd(s::Int64, N::Int64, k::Int64, norbs::Int64)
    function FaTerm(N, k)
        return (i) -> exp(-im * ((2 * pi * k * i) / N))
    end
    return checkstate(s, N, norbs, FaTerm(N, k))
end

"""
--------------------------------------------------------------------------------
representative functions (used in momentum-conservation)
Find the representative state for a given state. 
return 
-> The rep. state(as int) 
-> Number of translations needed to reach the rep. state from given state
-> Trace sign
--------------------------------------------------------------------------------
"""

function representative(a::Int64, abit::Vector{Int64}, N::Int64, norbs::Int64, updateRep::F) where {F<:Function}
    a_star = copy(a)
    l = 0
    tr_sign = 1
    rep = copy(a)

    for i in 1:N
        a_star = cyclebits(a_star, N, norbs)
        if a_star < rep
            rep, tr_sign, l = updateRep(a_star, abit, norbs, i)
        end
    end

    return rep, l, tr_sign
end

function representative_even(a::Int64, abit::Array{Int64,1}, N::Int64, norbs::Int64)
    function updateRepEven(a_star, abit, norbs, i)
        return a_star, (-1)^sum(abit[1:2*i*norbs]), i
    end
    return representative(a, abit, N, norbs, updateRepEven)
end

function representative_odd(a::Int64, abit::Vector{Int64}, N::Int64, norbs::Int64)
    function updateRepOdd(a_star, abit, norbs, i)
        return a_star, 1, i
    end
    return representative(a, abit, N, norbs, updateRepOdd)
end

"""
--------------------------------------------------------------------------------
formlist functions:
Form a list of all representative states in a given symmetry sector (nf/Sz/k/...)
return
-> List of rep. states
--------------------------------------------------------------------------------
"""

function formlist(N::Int64, norbs::Int64, condition::F) where {F<:Function}
    states = Int64[]
    endState = 4^(norbs * N) - 1
    for i in 0:endState
        condition(i) && push!(states, i)
    end
    return states
end

"""
--------------------------------------------------------------------------------
Conservation checking functions
--------------------------------------------------------------------------------
"""

function getFermionNumberConservation(nf::Int64)
    return (x) -> count_ones(x) == nf
end

function getSzConservation(Sz::Int64)
    return (x) -> countOddEvenDifference(x) == Sz
end

function getMomentumConservation(N::Int64, k::Int64, norbs::Int64, checkstateFunc::F) where {F<:Function}
    return (x) -> checkstateFunc(x, N, k, norbs) > 0
end

"""
--------------------------------------------------------------------------------
Functions to form list of states following given conservation laws
--------------------------------------------------------------------------------
"""

function formlist(N::Int64, norbs::Int64, condition::Symbol, conditionParameter::Vararg{T,M}) where {T,M}
    return formlist(Val(condition), N, norbs, conditionParameter...)
end

function formlist(::Val{T}, N, norbs, conditionParameter...) where {T}
    error(NotImplementedError())
end

function formlist(::Val{:nf}, N::Int64, norbs::Int64, nf::Int64)
    return formlist(N, norbs, getFermionNumberConservation(nf))
end

function formlist(::Val{:sz}, N::Int64, norbs::Int64, Sz::Int64)
    return formlist(N, norbs, getSzConservation(Sz))
end

function formlist(::Val{:nfsz}, N::Int64, norbs::Int64, nf::Int64, Sz::Int64)
    return formlist(N, norbs, (x) -> getFermionNumberConservation(nf)(x) && getSzConservation(Sz)(x))
end

function formlist(::Val{:nftr}, N::Int64, nf::Int64, k::Int64, norbs::Int64, checkstateFunc::F) where {F<:Function}
    basis = formlist(N, norbs, (x) -> getFermionNumberConservation(nf)(x) && getMomentumConservation(N, k, norbs, checkstateFunc)(x))
    return basis, checkstateFunc.(basis, N, k, norbs)
end

function formlist(T::Val{:nftr}, N::Int64, norbs::Int64, nf::Int64, k::Int64)
    if iseven(nf)
        return formlist(T, N, nf, k, norbs, checkstate_even)
    else
        return formlist(T, N, nf, k, norbs, checkstate_odd)
    end
end

function formlist(::Val{:nftrsz}, N::Int64, nf::Int64, Sz::Int64, k::Int64, norbs::Int64, checkstateFunc::F) where {F<:Function}
    basis = formlist(N, norbs,
        (x) -> getFermionNumberConservation(nf)(x) &&
                   getSzConservation(Sz)(x) &&
                   getMomentumConservation(N, k, norbs, checkstateFunc)(x)
    )
    return basis, checkstateFunc.(basis, N, k, norbs)
end

function formlist(T::Val{:nftrsz}, N::Int64, norbs::Int64, nf::Int64, Sz::Int64, k::Int64)
    if iseven(nf)
        return formlist(T, N, nf, Sz, k, norbs, checkstate_even)
    else
        return formlist(T, N, nf, Sz, k, norbs, checkstate_odd)
    end
end
module ED

using SparseArrays
using LinearAlgebra
using Arpack

include("Utils.jl")
include("Basis.jl")
export flip, cyclebits, findstate,
    checkstate, checkstate_even, checkstate_odd,
    representative, representative_even, representative_odd,
    formlist, getFermionNumberConservation, getSzConservation, getMomentumConservation
include("HubbardKanamori.jl")
export interactionTerms!, hoppingTerms!

end # module ED

@kwdef struct NotImplementedError <: Exception
    msg = "This method is currently not implemented"
end

function countOddEvenDifference(n::Int64)
    ODD_POSITION_ONES = 6148914691236517205
    EVEN_POSITION_ONES = -6148914691236517206
    return count_ones(n & ODD_POSITION_ONES) - count_ones(n & EVEN_POSITION_ONES)
end

"""
Internal method: Subject to change
"""
function _getBit(i, padding)
    return digits(i, base=2, pad=padding)
end
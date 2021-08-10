function hi!(x::Int)
    #this doesn't change input because Ints are immutable
    x = x + 2
    return x
end

function hi!(x::Array{T,1}) where T <: Real
    x .= x .+ 2
    return x
end
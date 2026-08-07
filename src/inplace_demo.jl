function hi!(x::Int)
    # This doesn't change the input because Ints are immutable.
    # This just makes the local name x refer to a new value.
    x = x + 2
    return x
end

function hi!(x::Array{T,1}) where T <: Real
    # .= changes the contents of the existing array.
    x .= x .+ 2
    return x
end

function liar_hi!(x::Array{T,1}) where T <: Real
    # Notice "=" versus ".="
    # "=" means "make the local name x refer to a new array"
    # ".=" means "change the things inside the existing array"
    x = x .+ 2
    return x
end

function sneaky_hi(x)
    # Exclamation marks are a hint to the developer.
    # They don't actually make Julia do anything special.
    x .= x .+ 2
    return x
end

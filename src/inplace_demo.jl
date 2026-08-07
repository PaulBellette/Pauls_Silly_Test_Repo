function hi!(x::Int)
    #this doesn't change input because Ints are immutable
    x = x + 2
    return x
end

function hi!(x::Array{T,1}) where T <: Real
    x .= x .+ 2
    return x
end

function liar_hi!(x::Array{T,1}) where T <: Real
    x = x .+ 2 #notice "=" versus ".=", "new thing called x in function" versus "set all the things in x to" 
    return x
end

function sneaky_hi(x)
    #exclaimation marks are a hint to developer, they don't actually do anything
    x .= x .+ 2
    return x
end

function sum_rows(data, n, m)
    s = zero(eltype(data))
    for i in 1:n
        for j in 1:m
            s += data[i,j]
        end
    end
    return s
end

function sum_columns(data, n, m)
    s = zero(eltype(data))
    for j in 1:m
        for i in 1:n
            s += data[i,j]
        end
    end
    return s
end

using Random

Random.seed!(0)
n, m = 500, 10^6
data = rand(n,m)

#dry run for precompilation
sum_columns(data, n, m)
sum_rows(data, n, m)

println("Summing with the column iteration moving fastest")
@time sum_columns(data, n, m)
println("Summing with the row iteration moving fastest")
@time sum_rows(data, n, m);

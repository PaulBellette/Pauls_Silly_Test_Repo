using Images, LinearAlgebra, Plots

function low_rank_approx(x::SVD{T}, k) where T <: Real
    #take a low rank (k) approx of the SVD of a matrix
    #via outer products
    return x.U[:, 1:k] * Diagonal(x.S[1:k]) * x.V[:, 1:k]'
end

image_path = joinpath(dirname(@__DIR__), "data/20210501_171236.jpg")

img_big = load(image_path)
img = Images.imresize(img_big, ratio = 1/4)

channels = channelview(img)
svd_r = svd(channels[1,:,:])
svd_g = svd(channels[2,:,:])
svd_b = svd(channels[3,:,:])

p1 = plot(svd_r.S, c = :red, label = "red singular values")
plot!(p1, svd_g.S, c = :green, label = "green singular values")
plot!(p1, svd_b.S, c = :blue, label = "blue singular values")

#recall trace norm identity

# ||A||f  = √sum(σ(A)^2)

#Frobenius norm => Pretend matrix is big vector and take l2 norm
# hints that truncated svd might contain a good representation on matrix

#in case you don't trust me...
norm(channels[1,:,:]) ≈ norm(svd_r.S)
norm(channels[2,:,:]) ≈ norm(svd_g.S)
norm(channels[3,:,:]) ≈ norm(svd_b.S)

approx_error = [norm(low_rank_approx(svd_r, i) - channels[1,:,:]) for i = 1 :10: length(svd_r.S)]

p2 = plot(approx_error./norm(channels[1,:,:]))

#what does it look like though...

reconstructed_images = []

for k = 1:10:60
    r_low = low_rank_approx(svd_r, k)
    g_low = low_rank_approx(svd_g, k)
    b_low = low_rank_approx(svd_b, k)
    push!(reconstructed_images, colorview(RGB, r_low, g_low, b_low))
end

output_plot = [plot(ri) for ri in reconstructed_images]

plot(output_plot...)

#So we can get a linear reduction in size by keeping the components of the 
#outer product U[:, 1:k], S[1:k], V[:, 1:k] and reconstructing via multiplication

#So...

@show sizeof(img_big)/1e6

#So to get approx 1% reconstruction error in Frobenius norm we need about half the singular values

percent_reduction_01 = findfirst(approx_error./norm(channels[1,:,:]) .< 0.01)/length(approx_error)

@show percent_reduction_01

#so the final "saved" image would be

@show (sizeof(img_big)/1e6) * percent_reduction_01

#not bad, but look at what jpeg gets...
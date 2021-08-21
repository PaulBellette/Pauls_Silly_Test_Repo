using SpecialFunctions, Plots, HCubature, LinearAlgebra, Random

function n_ball_volume_analytic(n_dims, r)
    V = π^(n_dims/2)*r^n_dims/gamma(n_dims/2 + 1)
    return V
end

function monte_carlo_sphere_test(n_samples, r, n_dims)
    #draw uniformly random points in n_space and see if distance to
    #origin is <= r
    Random.seed!(0)
    inside_ball = 0
    for i = 1 : n_samples
        if norm(rand(n_dims)) <= r
            inside_ball += 1
        end
    end
    ratio = inside_ball/n_samples
    volume = ratio*(2*r)^n_dims
    return volume
end

function cubature_sphere_volume(n_samples, r, n_dims)
    a = -ones(n_dims-1) .+ eps()
    b = ones(n_dims-1) .- eps()
    function f(x)
        inside_bit = r^2 - sum(x.^2)
        if inside_bit <= 0 
            return 0.0
        else
            return sqrt(r^2 - sum(x.^2))
        end
    end
    I, e = hcubature(f, a, b, maxevals = n_samples)
    return 2*I
end

r = 1.0
n_dims = 2:15

V = n_ball_volume_analytic.(n_dims, r)

p1 = plot(n_dims, V, xlabel = "dimension", ylabel = "volume", label = "analytic")

n_samples = Int(1e6)

V_mc = monte_carlo_sphere_test.(n_samples, r, n_dims)

plot!(p1, n_dims, V_mc, xlabel = "dimension", ylabel = "volume", label = "mc")

V_cubature = cubature_sphere_volume.(n_samples, r, n_dims)

plot!(p1, n_dims, V_cubature, ylim = [0, 6], xlabel = "dimension", ylabel = "volume", label = "cubature")

p2 = plot(n_dims, abs.(V_mc .- V)./V, yaxis = :log, xlabel = "dimension", ylabel = "relative abs error", label = "mc")

plot!(p2, n_dims, abs.(V_cubature .- V)./V, yaxis = :log, xlabel = "dimension", ylabel = "relative abs error", label = "cubature")

plot(p1, p2)
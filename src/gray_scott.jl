#!/usr/bin/env julia
# =============================================================================
# Gray–Scott reaction–diffusion → animated GIF
#
#   ∂u/∂t = Du·∇²u − u·v² + F·(1−u)
#   ∂v/∂t = Dv·∇²v + u·v² − (F+k)·v
#
# Two virtual chemicals U and V diffuse across a grid and react with each
# other. Depending on the feed rate F and kill rate k, the pattern that
# emerges from noise can look like coral, spots, worms, or spirals.
#
# Dependencies: NONE — pure Base Julia does the simulation. The only external
# tool used is the `ffmpeg` binary (must be on your PATH) to turn the raw
# frames into a looping GIF. If you don't have it: `brew install ffmpeg`,
# `sudo apt install ffmpeg`, or download from ffmpeg.org.
#
# Run with:  julia gray_scott.jl
# =============================================================================

using Random

# 5-point periodic Laplacian (in place)
function laplacian!(L::Matrix{Float64}, A::Matrix{Float64}, Ny::Int, Nx::Int)
    @inbounds for i in 1:Nx, j in 1:Ny
        jm = j == 1  ? Ny : j - 1
        jp = j == Ny ? 1  : j + 1
        im = i == 1  ? Nx : i - 1
        ip = i == Nx ? 1  : i + 1
        L[j, i] = A[jm, i] + A[jp, i] + A[j, im] + A[j, ip] - 4 * A[j, i]
    end
end

function palette(v::Float64)
    v = clamp(v, 0.0, 1.0)
    a = (0.5, 0.5, 0.5)
    b = (0.5, 0.5, 0.5)
    c = (1.0, 1.0, 1.0)
    d = (0.0, 0.10, 0.2)
    r = a[1] + b[1]*cos(2π*(c[1]*v + d[1]))
    g = a[2] + b[2]*cos(2π*(c[2]*v + d[2]))
    b_ = a[3] + b[3]*cos(2π*(c[3]*v + d[3]))
    return (round(UInt8, 255*clamp(r,0,1)), round(UInt8, 255*clamp(g,0,1)), round(UInt8, 255*clamp(b_,0,1)))
end

# write one frame as a binary PPM (no image packages needed)
function write_ppm(path::String, V::Matrix{Float64}, Ny::Int, Nx::Int)
    vmin, vmax = extrema(V)
    span = max(vmax - vmin, 1e-9)          # per-frame contrast stretch
    open(path, "w") do io
        write(io, "P6\n$Nx $Ny\n255\n")
        buf = Vector{UInt8}(undef, 3 * Nx)
        for j in 1:Ny
            idx = 1
            for i in 1:Nx
                r, g, b = palette((V[j, i] - vmin) / span)
                buf[idx] = r
                buf[idx+1] = g
                buf[idx+2] = b
                idx += 3
            end
            write(io, buf)
        end
    end
end

function main()
    Ny, Nx = 220, 220
    Du, Dv = 0.16, 0.08
    F, k = 0.0545, 0.0620   # "coral growth" — fills the grid, keeps a slow living texture, never fully dies out
    dt = 1.0
    steps = 9000
    save_every = 60                # -> 150 frames
    fps = 24
    gif_size = 480               # final (upscaled) pixel size
    max_colors = 128               # GIF palette size (lower = smaller file)

    outdir  = joinpath(@__DIR__, "frames")
    outfile = joinpath(@__DIR__, "reaction_diffusion.gif")
    isdir(outdir) && rm.(readdir(outdir, join=true))
    mkpath(outdir)

    # initial state: a handful of seed blobs
    Random.seed!(112358)
    U, V = ones(Ny, Nx), zeros(Ny, Nx)
    for _ in 1:4
        cy, cx, r = rand(40:Ny-40), rand(40:Nx-40), rand(15:25)
        @inbounds for i in 1:Nx, j in 1:Ny
            if (i - cx)^2 + (j - cy)^2 < r^2
                U[j, i] = 0.50
                V[j, i] = 0.25
            end
        end
    end

    Lu, Lv = similar(U), similar(V)
    Un, Vn = similar(U), similar(V)
    frame  = 0

    println("simulating $steps steps on a $(Ny)x$(Nx) grid...")
    for step in 1:steps
        laplacian!(Lu, U, Ny, Nx)
        laplacian!(Lv, V, Ny, Nx)
        @inbounds for i in 1:Nx, j in 1:Ny
            uvv = U[j, i] * V[j, i]^2
            Un[j, i] = U[j, i] + dt * (Du * Lu[j, i] - uvv + F * (1 - U[j, i]))
            Vn[j, i] = V[j, i] + dt * (Dv * Lv[j, i] + uvv - (F + k) * V[j, i])
        end
        U, Un = Un, U
        V, Vn = Vn, V

        if step % save_every == 0
            frame += 1
            write_ppm(joinpath(outdir, "frame_$(lpad(frame, 4, '0')).ppm"), V, Ny, Nx)
        end
    end
    println("done simulating — wrote $frame frames to $outdir")

    # assemble the GIF with ffmpeg (two-pass palette = clean colour)
    frames_glob = joinpath(outdir, "frame_%04d.ppm")
    palette_png = joinpath(outdir, "palette.png")
    vf_gen = "scale=$gif_size:$gif_size:flags=lanczos,palettegen=max_colors=$max_colors"
    vf_use = "[0:v]scale=$gif_size:$gif_size:flags=lanczos[x];[x][1:v]paletteuse=dither=none"

    try
        run(`ffmpeg -y -framerate $fps -i $frames_glob -vf $vf_gen $palette_png`)
        run(`ffmpeg -y -framerate $fps -i $frames_glob -i $palette_png -lavfi $vf_use -loop 0 $outfile`)
        println("done -> $outfile")
    catch e
        println()
        println("Couldn't run ffmpeg (frames are still saved in $outdir).")
        println("Install ffmpeg and be sure it's on your PATH, then run:")
        println("  ffmpeg -y -framerate $fps -i $frames_glob -vf \"$vf_gen\" $palette_png")
        println("  ffmpeg -y -framerate $fps -i $frames_glob -i $palette_png -lavfi \"$vf_use\" -loop 0 $outfile")
    end
end

main()

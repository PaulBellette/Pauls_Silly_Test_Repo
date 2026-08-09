#!/usr/bin/env julia
# =============================================================================
# Stable Fluids — colored ink swirling through a 2D velocity field
#
# Jos Stam's classic (1999) method: the velocity field is diffused, made
# divergence-free ("projected"), and advected through itself each step;
# three dye fields (R, G, B) ride along, diffusing and advecting through
# that velocity field, giving the look of colored ink stirred into water.
#
# Periodic (wrap-around) boundary, Jacobi relaxation for the linear solves —
# simpler to get right than reflecting walls / Gauss-Seidel, at the cost of
# a few more iterations to converge. This is the same structure used in
# Mike Ash's well-known "Fluid Simulation for Dummies" tutorial, adapted to
# wrap-around edges and three colour channels instead of one.
#
# Dependencies: NONE — pure Base Julia does the simulation. The only external
# tool used is the `ffmpeg` binary (must be on your PATH) to turn the raw
# frames into a looping GIF.
#
# Run with:  julia stable_fluids.jl
# =============================================================================

using Random

wrap(i::Int, N::Int) = mod(i - 1, N) + 1

# sum of the 4 periodic neighbours
function neighbor_sum(A, Ny, Nx)
    S = similar(A)
    @inbounds for i in 1:Nx, j in 1:Ny
        jm = j == 1  ? Ny : j - 1
        jp = j == Ny ? 1  : j + 1
        im = i == 1  ? Nx : i - 1
        ip = i == Nx ? 1  : i + 1
        S[j, i] = A[jm, i] + A[jp, i] + A[j, im] + A[j, ip]
    end
    return S
end

# implicit diffusion via Jacobi relaxation
function diffuse(x0, rate, dt, Ny, Nx, iters)
    a = dt * rate
    x = copy(x0)
    for _ in 1:iters
        S = neighbor_sum(x, Ny, Nx)
        @inbounds for i in 1:Nx, j in 1:Ny
            x[j, i] = (x0[j, i] + a * S[j, i]) / (1 + 4a)
        end
    end
    return x
end

# pressure projection: subtract the gradient to make velocity div-free
function project(vx, vy, Ny, Nx, iters)
    div = similar(vx)
    @inbounds for i in 1:Nx, j in 1:Ny
        im = i == 1  ? Nx : i - 1
        ip = i == Nx ? 1  : i + 1
        jm = j == 1  ? Ny : j - 1
        jp = j == Ny ? 1  : j + 1
        div[j, i] = -0.5 * ((vx[j, ip] - vx[j, im]) + (vy[jp, i] - vy[jm, i]))
    end
    p = zeros(Ny, Nx)
    for _ in 1:iters
        S = neighbor_sum(p, Ny, Nx)
        @inbounds for i in 1:Nx, j in 1:Ny
            p[j, i] = (div[j, i] + S[j, i]) / 4
        end
    end
    vx2, vy2 = copy(vx), copy(vy)
    @inbounds for i in 1:Nx, j in 1:Ny
        im = i == 1  ? Nx : i - 1
        ip = i == Nx ? 1  : i + 1
        jm = j == 1  ? Ny : j - 1
        jp = j == Ny ? 1  : j + 1
        vx2[j, i] -= 0.5 * (p[j, ip] - p[j, im])
        vy2[j, i] -= 0.5 * (p[jp, i] - p[jm, i])
    end
    return vx2, vy2
end

# semi-Lagrangian backward advection with bilinear interpolation
function advect(field, vx, vy, dt, Ny, Nx)
    out = similar(field)
    @inbounds for i in 1:Nx, j in 1:Ny
        px = i - dt * vx[j, i]
        py = j - dt * vy[j, i]
        bx = floor(Int, px); fx = px - bx
        by = floor(Int, py); fy = py - by
        i0 = wrap(bx, Nx);   i1 = wrap(bx + 1, Nx)
        j0 = wrap(by, Ny);   j1 = wrap(by + 1, Ny)
        out[j, i] = (1 - fx) * (1 - fy) * field[j0, i0] +
                    fx       * (1 - fy) * field[j0, i1] +
                    (1 - fx) * fy       * field[j1, i0] +
                    fx       * fy       * field[j1, i1]
    end
    return out
end

# add a Gaussian vortex (rotational velocity kick)
function add_vortex!(vx, vy, Ny, Nx, cy, cx, strength, radius)
    @inbounds for i in 1:Nx, j in 1:Ny
        dxp = i - cx
        dyp = j - cy
        falloff = exp(-(dxp^2 + dyp^2) / (2 * radius^2))
        vx[j, i] += -dyp * strength * falloff
        vy[j, i] +=  dxp * strength * falloff
    end
end

# add a Gaussian blob of dye
function add_blob!(field, Ny, Nx, cy, cx, radius, amp)
    @inbounds for i in 1:Nx, j in 1:Ny
        dxp = i - cx
        dyp = j - cy
        field[j, i] += amp * exp(-(dxp^2 + dyp^2) / (2 * radius^2))
    end
end

# write one RGB frame as a binary PPM
function write_ppm_rgb(path, R, G, B, Ny, Nx)
    open(path, "w") do io
        write(io, "P6\n$Nx $Ny\n255\n")
        buf = Vector{UInt8}(undef, 3 * Nx)
        for j in 1:Ny
            idx = 1
            for i in 1:Nx
                buf[idx]   = round(UInt8, 255 * clamp(R[j, i], 0.0, 1.0))
                buf[idx+1] = round(UInt8, 255 * clamp(G[j, i], 0.0, 1.0))
                buf[idx+2] = round(UInt8, 255 * clamp(B[j, i], 0.0, 1.0))
                idx += 3
            end
            write(io, buf)
        end
    end
end

function main()
    # parameters
    Ny, Nx = 150, 150
    visc, diffr = 0.015, 0.010 # velocity viscosity, dye diffusion
    dt = 0.2
    gs_iters = 25 # Jacobi relaxation sweeps per solve
    steps = 420
    save_every = 3 # -> 140 frames
    fps = 30
    gif_size = 480
    max_colors = 200
    dither = "bayer:bayer_scale=4"

    outdir = joinpath(@__DIR__, "frames")
    outfile = joinpath(@__DIR__, "stable_fluids.gif")
    isdir(outdir) && rm.(readdir(outdir, join=true))
    mkpath(outdir)

    vx, vy = zeros(Ny, Nx), zeros(Ny, Nx)
    dr, dg, db = zeros(Ny, Nx), zeros(Ny, Nx), zeros(Ny, Nx)

    # five anchor vortices in a cross pattern, alternating spin
    vortex_specs = [
        (Ny*0.30, Nx*0.30,  0.35, 18.0),
        (Ny*0.70, Nx*0.30, -0.35, 18.0),
        (Ny*0.30, Nx*0.70, -0.35, 18.0),
        (Ny*0.70, Nx*0.70,  0.35, 18.0),
        (Ny*0.50, Nx*0.50,  0.45, 22.0),
    ]
    for (cy, cx, s, r) in vortex_specs
        add_vortex!(vx, vy, Ny, Nx, cy, cx, s, r)
    end

    Random.seed!(3)
    for _ in 1:9
        cy, cx = 0.1*Ny + 0.8*Ny*rand(), 0.1*Nx + 0.8*Nx*rand()
        s = rand([-1, 1]) * (0.15 + 0.15*rand())
        r = 7 + 6*rand()
        add_vortex!(vx, vy, Ny, Nx, cy, cx, s, r)
    end

    # coloured ink blobs seeded near the vortices
    blob_specs = [
        (dr, Ny*0.30, Nx*0.30, 14.0, 1.0),
        (dg, Ny*0.70, Nx*0.30, 14.0, 1.0),
        (db, Ny*0.30, Nx*0.70, 14.0, 1.0),
        (dr, Ny*0.70, Nx*0.70, 14.0, 0.6),
        (dg, Ny*0.70, Nx*0.70, 14.0, 0.6),
        (db, Ny*0.50, Nx*0.50, 16.0, 1.0),
    ]
    for (field, cy, cx, r, amp) in blob_specs
        add_blob!(field, Ny, Nx, cy, cx, r, amp)
    end

    frame = 0
    println("simulating $steps steps on a $(Ny)x$(Nx) grid...")
    for s in 1:steps
        vx = diffuse(vx, visc, dt, Ny, Nx, gs_iters)
        vy = diffuse(vy, visc, dt, Ny, Nx, gs_iters)
        vx, vy = project(vx, vy, Ny, Nx, gs_iters)
        vx2 = advect(vx, vx, vy, dt, Ny, Nx)
        vy2 = advect(vy, vx, vy, dt, Ny, Nx)
        vx, vy = project(vx2, vy2, Ny, Nx, gs_iters)

        dr = advect(diffuse(dr, diffr, dt, Ny, Nx, gs_iters), vx, vy, dt, Ny, Nx)
        dg = advect(diffuse(dg, diffr, dt, Ny, Nx, gs_iters), vx, vy, dt, Ny, Nx)
        db = advect(diffuse(db, diffr, dt, Ny, Nx, gs_iters), vx, vy, dt, Ny, Nx)

        # a small random stir every so often, so it never fully settles
        if s % 55 == 0
            cy, cx = 0.15*Ny + 0.7*Ny*rand(), 0.15*Nx + 0.7*Nx*rand()
            sgn = rand([-1, 1]) * (0.15 + 0.13*rand())
            add_vortex!(vx, vy, Ny, Nx, cy, cx, sgn, 8 + 6*rand())
        end

        if s % save_every == 0
            frame += 1
            write_ppm_rgb(joinpath(outdir, "frame_$(lpad(frame, 4, '0')).ppm"), dr, dg, db, Ny, Nx)
        end
    end
    println("done simulating — wrote $frame frames to $outdir")

    # assemble the GIF with ffmpeg (two-pass palette + light dither)
    frames_glob = joinpath(outdir, "frame_%04d.ppm")
    palette_png = joinpath(outdir, "palette.png")
    vf_gen = "scale=$gif_size:$gif_size:flags=lanczos,palettegen=max_colors=$max_colors"
    vf_use = "[0:v]scale=$gif_size:$gif_size:flags=lanczos[x];[x][1:v]paletteuse=dither=$dither"

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

using InclusionAnalytics
using CairoMakie, MathTeXEngine, StaticArrays, Statistics, Printf

# Compare numerics vs analytics for either a circular or an elliptical inclusion
# with an out-of-plane strain rate ε̇zz activated.
let
    geometry = :elliptical   # :circular | :elliptical
    test     = :Duretz2026_7_oop
    nc       = 201
    params   = preset(test, geometry)
    f_anal   = geometry==:circular ? analytics_circle : analytics_ellipse

    # Numerics
    V, P, σ, X = Stokes2D(nc, test, params, f_anal)

    # Analytics
    Va, Pa, σa = analytics_stag(f_anal, X, params, geometry)

    # Errors
    Ve, Pe, σe = errors_stag(V, P, σ, Va, Pa, σa)

    # Coordinate normalisation
    Lx    = X.xv[end] - X.xv[1]
    scale = 2 / Lx
    xv, yv   = scale .* X.xv,  scale .* X.yv
    xc, yc   = scale .* X.xc,  scale .* X.yc
    xce, yce = scale .* X.xce, scale .* X.yce

    field_data = [
        (L"V_x",         xv,  yce, Va.x,   V.x,   Ve.x,   false),
        (L"V_y",         xce, yv,  Va.y,   V.y,   Ve.y,   false),
        (L"P",           xc,  yc,  Pa,     P,     Pe,     true ),
        (L"\sigma_{zz}", xc,  yc,  σa.zz,  σ.zz,  σe.zz,  true ),
    ]

    fs     = 26  
    fs_lbl = 34   

    set_theme!(theme_latexfonts())
    fig = Figure(size=(1300, 1500), fontsize=fs)
    Label(fig[0, 1:2], L"\mathrm{Analytical\ solution}", tellwidth=false, fontsize=fs_lbl)
    Label(fig[0, 3:4], L"\mathrm{Numerical\ solution}",  tellwidth=false, fontsize=fs_lbl)
    Label(fig[0, 5:6], L"\mathrm{Difference}",           tellwidth=false, fontsize=fs_lbl)

    for (row, (name, gx, gy, an, num, diff, share)) in enumerate(field_data)
        an_c  = an
        num_c = num

        mag = max(maximum(abs, an_c), maximum(abs, num_c), 1.0)
        crange = if share
            lo, hi  = extrema(vcat(vec(an_c), vec(num_c)))
            minspan = 0.02*mag
            if hi - lo < minspan
                c      = 0.5*(lo + hi)
                lo, hi = c - 0.5*minspan, c + 0.5*minspan
            end
            (colorrange = (lo, hi),)
        else
            NamedTuple()
        end

        dmax = max(maximum(abs, diff), 0.01*mag)

        bottom = (row == length(field_data))
        ax1 = Axis(fig[row, 1], aspect=DataAspect(), ylabel=L"y", xlabel=L"x",
                   xlabelsize=fs_lbl, ylabelsize=fs_lbl, xticklabelsize=20, yticklabelsize=20)
        ax2 = Axis(fig[row, 3], aspect=DataAspect(), xlabel=L"x",
                   xlabelsize=fs_lbl, xticklabelsize=20, yticklabelsize=20)
        ax3 = Axis(fig[row, 5], aspect=DataAspect(), xlabel=L"x",
                   xlabelsize=fs_lbl, xticklabelsize=20, yticklabelsize=20)

        hm1 = heatmap!(ax1, gx, gy, an_c;  colormap=(Reverse(:matter), 1), crange...)
        hm2 = heatmap!(ax2, gx, gy, num_c; colormap=(Reverse(:matter), 1), crange...)
        hm3 = heatmap!(ax3, gx, gy, diff;  colormap=(Reverse(:matter), 1), colorrange=(-dmax, dmax))

        Colorbar(fig[row, 2], hm1, label=name, labelsize=fs_lbl, ticklabelsize=20, width=14)
        Colorbar(fig[row, 4], hm2, label=name, labelsize=fs_lbl, ticklabelsize=20, width=14)
        Colorbar(fig[row, 6], hm3, width=14, ticklabelsize=20)

        hideydecorations!(ax2, grid=false)
        hideydecorations!(ax3, grid=false)
        if !bottom
            hidexdecorations!.((ax1, ax2, ax3), grid=false)
        end
    end

    colgap!(fig.layout, 10)
    rowgap!(fig.layout, 10)
    display(fig)

    # Save figure
    fname = "./figures/comparison_$(geometry)_nc$(nc)_$(test).png"
    save(fname, fig, px_per_unit=2)
    println("Saved: $fname")

    τa_xx, τ_xx = σa.xx .+ Pa, σ.xx .+ P
    τa_yy, τ_yy = σa.yy .+ Pa, σ.yy .+ P
    τa_zz, τ_zz = σa.zz .+ Pa, σ.zz .+ P
    τa_xy, τ_xy = σa.xy, σ.xy

    # Second invariant of the deviatoric stress (τxy averaged from vertices to cell centres)
    τaII = sqrt.(0.5*(τa_xx.^2 .+ τa_yy.^2 .+ τa_zz.^2) .+
                 (0.25*(τa_xy[1:end-1,1:end-1] .+ τa_xy[1:end-1,2:end] .+ τa_xy[2:end,1:end-1] .+ τa_xy[2:end,2:end])).^2)
    τII  = sqrt.(0.5*(τ_xx.^2  .+ τ_yy.^2  .+ τ_zz.^2 ) .+
                 (0.25*(τ_xy[1:end-1,1:end-1]  .+ τ_xy[1:end-1,2:end]  .+ τ_xy[2:end,1:end-1]  .+ τ_xy[2:end,2:end] )).^2)

    devstress_data = [
        (L"\tau_{xx}",  xc, yc, τa_xx, τ_xx, abs.(τ_xx .- τa_xx), true),
        (L"\tau_{yy}",  xc, yc, τa_yy, τ_yy, abs.(τ_yy .- τa_yy), true),
        (L"\tau_{xy}",  xv, yv, τa_xy, τ_xy, abs.(τ_xy .- τa_xy), true),
        (L"\tau_{zz}",  xc, yc, τa_zz, τ_zz, abs.(τ_zz .- τa_zz), true),
        (L"\tau_{II}",  xc, yc, τaII,  τII,  abs.(τII .- τaII),   true),
    ]

    fig2 = Figure(size=(1300, 1850), fontsize=fs)
    Label(fig2[0, 1:2], L"\mathrm{Analytical\ solution}", tellwidth=false, fontsize=fs_lbl)
    Label(fig2[0, 3:4], L"\mathrm{Numerical\ solution}",  tellwidth=false, fontsize=fs_lbl)
    Label(fig2[0, 5:6], L"\mathrm{Difference}",           tellwidth=false, fontsize=fs_lbl)

    for (row, (name, gx, gy, an, num, diff, share)) in enumerate(devstress_data)
        an_c  = an
        num_c = num

        mag = max(maximum(abs, an_c), maximum(abs, num_c), 1.0)
        crange = if share
            lo, hi  = extrema(vcat(vec(an_c), vec(num_c)))
            minspan = 0.02*mag
            if hi - lo < minspan
                c      = 0.5*(lo + hi)
                lo, hi = c - 0.5*minspan, c + 0.5*minspan
            end
            (colorrange = (lo, hi),)
        else
            NamedTuple()
        end

        dmax = max(maximum(abs, diff), 0.01*mag)

        bottom = (row == length(devstress_data))
        ax1 = Axis(fig2[row, 1], aspect=DataAspect(), ylabel=L"y", xlabel=L"x",
                   xlabelsize=fs_lbl, ylabelsize=fs_lbl, xticklabelsize=20, yticklabelsize=20)
        ax2 = Axis(fig2[row, 3], aspect=DataAspect(), xlabel=L"x",
                   xlabelsize=fs_lbl, xticklabelsize=20, yticklabelsize=20)
        ax3 = Axis(fig2[row, 5], aspect=DataAspect(), xlabel=L"x",
                   xlabelsize=fs_lbl, xticklabelsize=20, yticklabelsize=20)

        hm1 = heatmap!(ax1, gx, gy, an_c;  colormap=(Reverse(:matter), 1), crange...)
        hm2 = heatmap!(ax2, gx, gy, num_c; colormap=(Reverse(:matter), 1), crange...)
        hm3 = heatmap!(ax3, gx, gy, diff;  colormap=(Reverse(:matter), 1), colorrange=(-dmax, dmax))

        Colorbar(fig2[row, 2], hm1, label=name, labelsize=fs_lbl, ticklabelsize=20, width=14)
        Colorbar(fig2[row, 4], hm2, label=name, labelsize=fs_lbl, ticklabelsize=20, width=14)
        Colorbar(fig2[row, 6], hm3, width=14, ticklabelsize=20)

        hideydecorations!(ax2, grid=false)
        hideydecorations!(ax3, grid=false)
        if !bottom
            hidexdecorations!.((ax1, ax2, ax3), grid=false)
        end
    end

    colgap!(fig2.layout, 10)
    rowgap!(fig2.layout, 10)
    display(fig2)

    fname2 = "./figures/deviatoric_stresses_$(geometry)_nc$(nc)_$(test).png"
    save(fname2, fig2, px_per_unit=2)
    println("Saved: $fname2")
end

using InclusionAnalytics
using CairoMakie, MathTeXEngine, StaticArrays, Statistics, Printf

# Compare numerics vs analytics for the circular inclusion with an
# out-of-plane strain rate ε̇zz activated (only analytics_circle supports it).
let
    geometry = :circular
    test     = :Duretz2026_7_oop
    nc       = 301           
    params   = preset(test, geometry)
    f_anal   = analytics_circle

    # Numerics
    V, P, σ, X = Stokes2D(nc, test, params, f_anal)

    # Analytics
    Va, Pa, σa = analytics_stag(f_anal, X, params, geometry)

    # Errors
    Ve, Pe, σe = errors_stag(V, P, σ, Va, Pa, σa)

    # Coordinate normalisation: map the box onto [-1, 1] (physics untouched)
    Lx    = X.xv[end] - X.xv[1]
    scale = 2 / Lx
    xv, yv   = scale .* X.xv,  scale .* X.yv
    xc, yc   = scale .* X.xc,  scale .* X.yc
    xce, yce = scale .* X.xce, scale .* X.yce

    # Visualisation (velocities + pressure + out-of-plane stress)
    #   last entry: share a common colour range between analytical & numerical
    field_data = [
        (L"V_x",         xv,  yce, Va.x,   V.x,   Ve.x,   false),
        (L"V_y",         xce, yv,  Va.y,   V.y,   Ve.y,   false),
        (L"P",           xc,  yc,  Pa,     P,     Pe,     true ),
        (L"\sigma_{zz}", xc,  yc,  σa.zz,  σ.zz,  σe.zz,  true ),
    ]

    fs     = 26   # base font size
    fs_lbl = 34   # column titles / axis labels

    set_theme!(theme_latexfonts())
    fig = Figure(size=(1300, 1500), fontsize=fs)
    Label(fig[0, 1:2], L"\mathrm{Analytical\ solution}", tellwidth=false, fontsize=fs_lbl)
    Label(fig[0, 3:4], L"\mathrm{Numerical\ solution}",  tellwidth=false, fontsize=fs_lbl)
    Label(fig[0, 5:6], L"\mathrm{Difference}",           tellwidth=false, fontsize=fs_lbl)

    for (row, (name, gx, gy, an, num, diff, share)) in enumerate(field_data)
        an_c  = an
        num_c = num

        # common colour range for the analytical/numerical columns, widened to at
        # least 2% of the field's own magnitude so a (near-)uniform field doesn't
        # collapse into Float32-precision-limited patches
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

        # diff panel: zero-centred, widened the same way (relative to the field's
        # own magnitude, not the-often-tiny-diff itself) so it stays readable
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
end

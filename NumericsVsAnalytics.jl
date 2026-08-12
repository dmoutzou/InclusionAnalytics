using InclusionAnalytics
using CairoMakie, MathTeXEngine, StaticArrays, Statistics, Printf

# Compare numerics vs analytics for either a circular or an elliptical inclusion.
let
    geometry = :elliptical   # :circular | :elliptical
    test     = :Duretz2026_5_mixed   # :Schmid2003 | :Duretz2026_1_ps | :Duretz2026_2_ss | :Duretz2026_3_exp | :Duretz2026_4_comp | :Duretz2026_5_mixed | :Duretz2026_6_angled | :Duretz2026_7_oop
    nc       = 201           # resolution
    params   = preset(test, geometry)
    f_anal   = geometry==:circular ? analytics_circle : analytics_ellipse

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

    # Visualisation (velocities + pressure only)
    #   last entry: share a common colour range between analytical & numerical
    field_data = [
        (L"V_x", xv,  yce, Va.x, V.x, Ve.x, false),
        (L"V_y", xce, yv,  Va.y, V.y, Ve.y, false),
        (L"P",   xc,  yc,  Pa,   P,   Pe  , true ),
    ]

    fs     = 26   # base font size
    fs_lbl = 34   # column titles / axis labels

    set_theme!(theme_latexfonts())
    fig = Figure(size=(1300, 1150), fontsize=fs)
    Label(fig[0, 1:2], L"\mathrm{Analytical\ solution}", tellwidth=false, fontsize=fs_lbl)
    Label(fig[0, 3:4], L"\mathrm{Numerical\ solution}",  tellwidth=false, fontsize=fs_lbl)
    Label(fig[0, 5:6], L"\mathrm{Difference}",           tellwidth=false, fontsize=fs_lbl)

    for (row, (name, gx, gy, an, num, diff, share)) in enumerate(field_data)
        an_c  = an  
        num_c = num 

        # common colour range for the analytical/numerical columns
        crange = share ? (colorrange = extrema(vcat(vec(an_c), vec(num_c))),) : NamedTuple()

        bottom = (row == length(field_data))
        ax1 = Axis(fig[row, 1], aspect=DataAspect(), ylabel=L"y", xlabel=L"x",
                   xlabelsize=fs_lbl, ylabelsize=fs_lbl, xticklabelsize=20, yticklabelsize=20)
        ax2 = Axis(fig[row, 3], aspect=DataAspect(), xlabel=L"x",
                   xlabelsize=fs_lbl, xticklabelsize=20, yticklabelsize=20)
        ax3 = Axis(fig[row, 5], aspect=DataAspect(), xlabel=L"x",
                   xlabelsize=fs_lbl, xticklabelsize=20, yticklabelsize=20)

        hm1 = heatmap!(ax1, gx, gy, an_c;  colormap=(Reverse(:matter), 1), crange...)
        hm2 = heatmap!(ax2, gx, gy, num_c; colormap=(Reverse(:matter), 1), crange...)
        hm3 = heatmap!(ax3, gx, gy, diff;  colormap=(Reverse(:matter), 1))

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
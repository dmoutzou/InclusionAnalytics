using InclusionAnalytics
using CairoMakie, MathTeXEngine, StaticArrays, Statistics

set_theme!(theme_latexfonts())

let
    geometry1 = :circular
    geometry2 = :elliptical
    test      = :Duretz2026_5_mixed # :Schmid2003 | :Duretz2026_1_ps | :Duretz2026_2_ss | :Duretz2026_3_exp | :Duretz2026_4_comp | :Duretz2026_5_mixed | :Duretz2026_6_angled
    nc        = 501
    params1 = preset(test, geometry1)
    params2 = preset(test, geometry2)

    f_anal1 = analytics_circle
    f_anal2 = analytics_ellipse

    Lx, Ly   = 1.0, 1.0
    Δx, Δy   = Lx/nc, Ly/nc
    xce, yce = LinRange(-Lx/2-Δx/2, Lx/2+Δx/2, nc+2), LinRange(-Ly/2-Δy/2, Ly/2+Δy/2, nc+2)
    xc, yc   = LinRange(-Lx/2+Δx/2, Lx/2-Δx/2, nc),   LinRange(-Ly/2+Δy/2, Ly/2-Δy/2, nc)
    xv, yv   = LinRange(-Lx/2, Lx/2, nc+1),           LinRange(-Ly/2, Ly/2, nc+1)
    X  = (xv=xv, yv=yv, xce=xce, yce=yce, xc=xc, yc=yc)

    Va1, Pa1, σa1 = analytics_stag(f_anal1, X, params1, geometry1)
    Va2, Pa2, σa2 = analytics_stag(f_anal2, X, params2, geometry2)

    # Purely cosmetic rescale for plotting: stretches axis labels to [-1,1]
    # without touching the physics or the inclusion-to-domain ratio.
    scale = 2 / Lx   # assumes Lx == Ly; use separate scales if not
    xc_p, yc_p   = scale .* xc,  scale .* yc
    xv_p, yv_p   = scale .* xv,  scale .* yv
    xce_p, yce_p = scale .* xce, scale .* yce

    row1 = [
        (L"P",    xc_p,  yc_p,  Pa1   ),
        (L"\sigma_{xx}", xc_p,  yc_p,  σa1.xx),
        (L"\sigma_{yy}", xc_p,  yc_p,  σa1.yy),
        (L"\sigma_{xy}", xv_p,  yv_p,  σa1.xy),
    ]

    row2 = [
        (L"P",    xc_p,  yc_p,  Pa2   ),
        (L"\sigma_{xx}", xc_p,  yc_p,  σa2.xx),
        (L"\sigma_{yy}", xc_p,  yc_p,  σa2.yy),
        (L"\sigma_{xy}", xv_p,  yv_p,  σa2.xy),
    ]

    fig = Figure(size=(1200, 700), fontsize=12)

    for (col, (name, gx, gy, field)) in enumerate(row1)
        field_c = field .- mean(field)
        ax = Axis(fig[1, 2col-1], aspect=DataAspect(), title=name, titlesize=25)
        hm = heatmap!(ax, gx, gy, field_c; colormap=(Reverse(:matter), 1))
        Colorbar(fig[1, 2col], hm, width=12)
    end

    for (col, (name, gx, gy, field)) in enumerate(row2)
        field_c = field .- mean(field)
        ax = Axis(fig[2, 2col-1], aspect=DataAspect())
        hm = heatmap!(ax, gx, gy, field_c; colormap=(Reverse(:matter), 1))
        Colorbar(fig[2, 2col], hm, width=12)
    end

    Label(fig[1, 0], "Circular",   rotation=pi/2, tellheight=false, fontsize=20)
    Label(fig[2, 0], "Elliptical", rotation=pi/2, tellheight=false, fontsize=20)

    display(fig)

    fname = "./figures/analytics_circle_vs_ellipse_nc$(nc)_$(test).png"
    save(fname, fig, px_per_unit=2)
    println("Saved: $fname")
end
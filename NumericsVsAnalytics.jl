using InclusionAnalytics
using CairoMakie, MathTeXEngine, StaticArrays, Statistics, Printf

# Compare numerics vs analytics for either a circular or an elliptical inclusion.
let
    geometry = :circular   # :circular | :elliptical
    test     = :Duretz2026_5_mixed   # :Schmid2003 | :Duretz2026_1_ps | :Duretz2026_2_ss | :Duretz2026_3_exp | :Duretz2026_4_comp | :Duretz2026_5_mixed | :Duretz2026_6_angled
    nc       = 301           # resolution
    params   = preset(test, geometry)
    f_anal   = geometry==:circular ? analytics_circle : analytics_ellipse
    
    # Numerics
    V, P, σ, X = Stokes2D(nc, test, params, f_anal)
    
    # Analytics
    Va, Pa, σa = analytics_stag(f_anal, X, params, geometry)

    # Errors
    Ve, Pe, σe = errors_stag(V, P, σ, Va, Pa, σa)

    # Visualisation
    field_data = [
        ("Vx",  X.xv,  X.yce, Va.x,  V.x,  Ve.x ),
        ("Vy",  X.xce, X.yv,  Va.y,  V.y,  Ve.y ),
        ("P",   X.xc,  X.yc,  Pa,    P,    Pe   ),
        ("Sxx", X.xc,  X.yc,  σa.xx, σ.xx, σe.xx),
        ("Syy", X.xc,  X.yc,  σa.yy, σ.yy, σe.yy),
        ("Sxy", X.xv,  X.yv,  σa.xy, σ.xy, σe.xy),
    ]

    fig = Figure(size=(1100, 1700), fontsize=12)
    Label(fig[-1, 1:6], "geometry=$geometry  |  test=$test  |  nc=$nc",
          tellwidth=false, fontsize=15, font=:bold)
    Label(fig[0, 1:2], "Analytical", tellwidth=false, fontsize=14, font=:bold)
    Label(fig[0, 3:4], "Numerical",  tellwidth=false, fontsize=14, font=:bold)
    Label(fig[0, 5:6], "Difference", tellwidth=false, fontsize=14, font=:bold)

    for (row, (name, gx, gy, an, num, diff)) in enumerate(field_data)
        an_c  = an  .- mean(an)
        num_c = num .- mean(num)

        ax1 = Axis(fig[row, 1], aspect=DataAspect(), ylabel=name)
        hm1 = heatmap!(ax1, gx, gy, an_c;  colormap=(Reverse(:matter), 1))
        Colorbar(fig[row, 2], hm1, width=12)

        ax2 = Axis(fig[row, 3], aspect=DataAspect())
        hm2 = heatmap!(ax2, gx, gy, num_c; colormap=(Reverse(:matter), 1))
        Colorbar(fig[row, 4], hm2, width=12)

        ax3 = Axis(fig[row, 5], aspect=DataAspect())
        hm3 = heatmap!(ax3, gx, gy, diff;  colormap=(Reverse(:matter), 1))
        Colorbar(fig[row, 6], hm3, width=12)
    end
    display(fig)

    # Save figure
    fname = "./figures/comparison_$(geometry)_nc$(nc)_$(test).png"
    save(fname, fig, px_per_unit=2)
    println("Saved: $fname")
end

using InclusionAnalytics
using CairoMakie, MathTeXEngine, StaticArrays, Statistics, Printf

# Compare numerics vs analytics for either a circular or an elliptical inclusion.
# Stokes2D (src/Numerics_v2.jl) only solves the numerical problem; the analytical
# solution is evaluated here, on the same grids.
let
    geometry = :circular   # :circular | :elliptical
    test     = :Schmid2003   # :Schmid2003 | :Duretz2026_1 | :Duretz2026_2 | :Duretz2026_3
    nc       = 101           # resolution
    a_target = 0.4           # only used when geometry == :elliptical (normalised semi-major axis)

    V, Pt, S, X, params, scale = Stokes2D(nc, test; geometry=geometry, a_target=a_target)

    Va, Pa, τa = eval_analytics_stag(f_anal, X, params, geometry, scale)

    # Errors
    dVx = abs.(V.x .- Va.x);  dVy = abs.(V.y .- Va.y);  dP = abs.(Pt .- Pa)
    dSxx = abs.(S.xx .- τa.xx);  dSyy = abs.(S.yy .- τa.yy)
    dSzz = abs.(S.zz .- τa.zz);  dSxy = abs.(S.xy .- τa.xy)
    @printf("mean|dVx| = %1.3e   mean|dVy| = %1.3e   mean|dP| = %1.3e\n", mean(dVx), mean(dVy), mean(dP))
    @printf("mean|dSxx| = %1.3e   mean|dSyy| = %1.3e   mean|dSzz| = %1.3e   mean|dSxy| = %1.3e\n",
            mean(dSxx), mean(dSyy), mean(dSzz), mean(dSxy))

    field_data = [
        ("Vx",  X.xv,  X.yce, Va.x,  V.x,  dVx ),
        ("Vy",  X.xce, X.yv,  Va.y,  V.y,  dVy ),
        ("P",   X.xc,  X.yc,  Pa,    Pt,   dP  ),
        ("Sxx", X.xc,  X.yc,  τa.xx, S.xx, dSxx),
        ("Syy", X.xc,  X.yc,  τa.yy, S.yy, dSyy),
        ("Sxy", X.xv,  X.yv,  τa.xy, S.xy, dSxy),
    ]

    fig = Figure(size=(1100, 1700), fontsize=12)
    Label(fig[-1, 1:6], "geometry=$geometry  |  test=$test  |  nc=$nc" *
          (geometry == :elliptical ? "  |  a_target=$a_target" : ""),
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
        hm2 = heatmap!(ax2, gx, gy, num_c; colormap=(Reverse(:matter), 1), colorrange=extrema(an_c))
        Colorbar(fig[row, 4], hm2, width=12)

        ax3 = Axis(fig[row, 5], aspect=DataAspect())
        hm3 = heatmap!(ax3, gx, gy, diff;  colormap=(Reverse(:matter), 1))
        Colorbar(fig[row, 6], hm3, width=12)
    end

    fname = "./figures/comparison_$(geometry)_nc$(nc)_$(test).png"
    save(fname, fig, px_per_unit=2)
    display(fig)
    println("Saved: $fname")
end

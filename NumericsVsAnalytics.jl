using InclusionAnalytics
using CairoMakie, MathTeXEngine, StaticArrays, Statistics, Printf

# Compare numerics vs analytics for either a circular or an elliptical inclusion.
# Stokes2D (src/Numerics_v2.jl) only solves the numerical problem; the analytical
# solution is evaluated here, on the same grids.
let
    geometry = :elliptical   # :circular | :elliptical
    test     = :Schmid2003   # :Schmid2003 | :Duretz2026_1 | :Duretz2026_2 | :Duretz2026_3
    nc       = 301           # resolution
    a_target = 0.4           # only used when geometry == :elliptical (normalised semi-major axis)

    V, Pt, S, X, params, scale = Stokes2D(nc, test; geometry=geometry, a_target=a_target)

    # Analytical reference on the same staggered grids as the numerics
    Va = (
        x = zeros(size(V.x)),
        y = zeros(size(V.y)),
    )
    Pa = zeros(size(Pt))
    Sa = (
        xx = zeros(size(Pt)),
        yy = zeros(size(Pt)),
        zz = zeros(size(Pt)),
        xy = zeros(size(S.xy)),
    )
    for I in CartesianIndices(V.x)
        i, j = I[1], I[2]
        Va.x[i,j] = f_anal(@SVector([X.xv[i]; X.yce[j]]); params=params, geometry=geometry, S=scale).V[1]
    end
    for I in CartesianIndices(V.y)
        i, j = I[1], I[2]
        Va.y[i,j] = f_anal(@SVector([X.xce[i]; X.yv[j]]); params=params, geometry=geometry, S=scale).V[2]
    end
    for I in CartesianIndices(Pt)
        i, j = I[1], I[2]
        sol = f_anal(@SVector([X.xc[i]; X.yc[j]]); params=params, geometry=geometry, S=scale)
        Pa[i,j]    = sol.p
        Sa.xx[i,j] = sol.τ[1,1] - sol.p
        Sa.yy[i,j] = sol.τ[2,2] - sol.p
        Sa.zz[i,j] = sol.τzz    - sol.p
    end
    for I in CartesianIndices(S.xy)
        i, j = I[1], I[2]
        Sa.xy[i,j] = f_anal(@SVector([X.xv[i]; X.yv[j]]); params=params, geometry=geometry, S=scale).τ[1,2]
    end

    # Errors
    dVx = abs.(V.x .- Va.x);  dVy = abs.(V.y .- Va.y);  dP = abs.(Pt .- Pa)
    dSxx = abs.(S.xx .- Sa.xx);  dSyy = abs.(S.yy .- Sa.yy)
    dSzz = abs.(S.zz .- Sa.zz);  dSxy = abs.(S.xy .- Sa.xy)
    @printf("mean|dVx| = %1.3e   mean|dVy| = %1.3e   mean|dP| = %1.3e\n", mean(dVx), mean(dVy), mean(dP))
    @printf("mean|dSxx| = %1.3e   mean|dSyy| = %1.3e   mean|dSzz| = %1.3e   mean|dSxy| = %1.3e\n",
            mean(dSxx), mean(dSyy), mean(dSzz), mean(dSxy))

    field_data = [
        ("Vx",  X.xv,  X.yce, Va.x,  V.x,  dVx ),
        ("Vy",  X.xce, X.yv,  Va.y,  V.y,  dVy ),
        ("P",   X.xc,  X.yc,  Pa,    Pt,   dP  ),
        ("Sxx", X.xc,  X.yc,  Sa.xx, S.xx, dSxx),
        ("Syy", X.xc,  X.yc,  Sa.yy, S.yy, dSyy),
        ("Sxy", X.xv,  X.yv,  Sa.xy, S.xy, dSxy),
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
        hm2 = heatmap!(ax2, gx, gy, num_c; colormap=(Reverse(:matter), 1))
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

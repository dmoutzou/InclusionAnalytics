# Include packages and scripts
using ExactFieldSolutions
include("../src/InclusionAnalytics.jl")


define_params(params) = (
    ηm=params.ηm, ηi=params.ηi, ξm=params.ξm, ξi=params.ξi,
    ri=params.r2, t=params.t, α=params.α,
    ε̇=params.ε̇, γ̇=params.γ̇, ζ̇=params.ζ̇, ε̇zz=params.ε̇zz,
)
analytics_circle(X; params)  = Stokes2D_Moutzouris_circle(X;  params=define_params(params))
analytics_ellipse(X; params) = Stokes2D_Moutzouris_ellipse(X; params=define_params(params))

set_theme!(theme_latexfonts())

function make_grid(nc; Lx=1.0, Ly=1.0)
    Δx, Δy = Lx/nc, Ly/nc
    xce = LinRange(-Lx/2-Δx/2, Lx/2+Δx/2, nc+2); yce = LinRange(-Ly/2-Δy/2, Ly/2+Δy/2, nc+2)
    xc  = LinRange(-Lx/2+Δx/2, Lx/2-Δx/2, nc);   yc  = LinRange(-Ly/2+Δy/2, Ly/2-Δy/2, nc)
    xv  = LinRange(-Lx/2, Lx/2, nc+1);           yv  = LinRange(-Ly/2, Ly/2, nc+1)
    return (xv=xv, yv=yv, xce=xce, yce=yce, xc=xc, yc=yc)
end

let
    geometry = :elliptical
    f_anal   = analytics_ellipse
    nc       = 501

    # single case with mixed loading
    test     = :Moutzouris2026_5_mixed
    tvalues  = [1.01, 2.0, 3.0, 4.0]

    # physical inclusion radius: hard-coded in preset() for every Duretz case (src/PhysicalParameters.jl)
    ri_phys  = 0.1

    # fixed unit box, as in src/Numerics.jl: r1/sc, r2/sc are the ellipse's
    Lx, Ly = 1.0, 1.0
    X      = make_grid(nc; Lx=Lx, Ly=Ly)
    scale  = 2 / Lx
    xc_p, yc_p = scale .* X.xc, scale .* X.yc

    fig   = Figure(size = (500 * length(tvalues), 700), fontsize = 18)
    cmap  = (Reverse(:matter), 1)
    ncols = 2 * length(tvalues)   # [axis | colorbar] per panel

    # constant parameters (t, r1, r2, sc vary across panels; the rest are fixed)
    params0 = preset(test, geometry)
    ηm, ηi, ξm, ξi, ri, r2, t0, α, sc0, ε̇, γ̇, ζ̇, ε̇zz = params0

    Label(fig[1, 1:ncols],
          L"\mu_h=%$(fmtnum(ηm))\;\;\mu_i=%$(fmtnum(ηi))\;\;\xi_h=%$(fmtnum(ξm))\;\;\xi_i=%$(fmtnum(ξi))\;\;\varepsilon=%$(fmtnum(ε̇))\;\;\gamma=%$(fmtnum(γ̇))\;\;\zeta=%$(fmtnum(ζ̇))\;\;\alpha=%$(fmtnum(α))",
          fontsize = 26, color = :gray25)

    panels = map(tvalues) do t
        params = preset(test, geometry)
        r1, r2 = ellipse_axes(t)
        sc     = r2 / ri_phys
        params = merge(params, (t = t, r1 = r1/sc, r2 = r2/sc, sc = sc))
        _, Pa, _ = analytics_stag(f_anal, X, params, geometry)
        return Pa
    end
    crange = (minimum(minimum, panels), maximum(maximum, panels))

    for (idx, t) in enumerate(tvalues)
        base = (idx - 1) * 2             # per panel: [axis | colorbar]
        Pa   = panels[idx]

        ax = Axis(fig[2, base+1], aspect = DataAspect(),
                  title = L"t=%$(fmtnum(t))", titlesize = 30, titlegap = 5,
                  xlabel = L"x", ylabel = L"y", xlabelsize = 28, ylabelsize = 28,
                  xticklabelsize = 18, yticklabelsize = 18)
        hm = heatmap!(ax, xc_p, yc_p, Pa; colorrange = crange, colormap = cmap)

        if idx == length(tvalues)
            Colorbar(fig[2, base+2], hm, width = 12, label = L"p", labelsize = 34, ticklabelsize = 18)
        else
            Colorbar(fig[2, base+2], hm, width = 12, ticklabelsize = 18)
        end
    end

    colgap!(fig.layout, 8)
    rowgap!(fig.layout, 1, 4)   

    #display(fig)

    mkpath("figures")
    fname = "./figures/Fig06_Analytics_mixed_aspect_ratio.png"
    save(fname, fig, px_per_unit = 2)
    println("Saved: $fname")
    fig
end
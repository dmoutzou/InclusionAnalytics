using InclusionAnalytics
using CairoMakie, MathTeXEngine, StaticArrays, Statistics
using ExactFieldSolutions

set_theme!(theme_latexfonts())

# Schmid & Podladchikov (2003) incompressible elliptical-inclusion solution,
# evaluated with the same interface as analytics_ellipse for direct comparison
function analytics_schmid_ellipse(X; params)
    ηm, ηi, ξm, ξi, ri, r2, t, α, sc, ε̇, γ̇, ζ̇, ε̇zz = params
    mc  = ηi / ηm
    Ζ   = InclusionAnalytics.to_zeta(sc .* X)
    sol = Analytics_Schmid2003(Ζ; ηm=ηm, mc=mc, γ̇=γ̇, ε̇=ε̇, α=α, t=t)
    return (V   = @SVector([sol.V[1]/sc, sol.V[2]/sc]),
            p   = sol.p,
            τ   = sol.τ,
            τzz = 0.0)
end

function make_grid(nc; Lx=1.0, Ly=1.0)
    Δx, Δy = Lx/nc, Ly/nc
    xce = LinRange(-Lx/2-Δx/2, Lx/2+Δx/2, nc+2); yce = LinRange(-Ly/2-Δy/2, Ly/2+Δy/2, nc+2)
    xc  = LinRange(-Lx/2+Δx/2, Lx/2-Δx/2, nc);   yc  = LinRange(-Ly/2+Δy/2, Ly/2-Δy/2, nc)
    xv  = LinRange(-Lx/2, Lx/2, nc+1);           yv  = LinRange(-Ly/2, Ly/2, nc+1)
    return (xv=xv, yv=yv, xce=xce, yce=yce, xc=xc, yc=yc)
end

let
    geometry = :elliptical                 # switch to :circular if you prefer
    f_anal   = geometry == :circular ? analytics_circle : analytics_ellipse
    nc       = 501

    Lx, Ly = 1.0, 1.0
    X      = make_grid(nc; Lx=Lx, Ly=Ly)

    scale = 2 / Lx                          # cosmetic axis rescale to [-1, 1]
    xc_p, yc_p = scale .* X.xc, scale .* X.yc

    base   = preset(:Duretz2026_1_ps, geometry)
    ξvals  = [0.1, 1.0e0, 1.0e3, 1.0e8]     # matrix (= inclusion) bulk modulus, last ~ incompressible
    incomp = 1.0e8                          # threshold to label as ->infty

    coltitle(ξ) = ξ ≥ incomp ? L"\xi_h=\xi_i\to\infty" : L"\xi_h=\xi_i=%$(fmtnum(ξ))"

    fig  = Figure(size = (1850, 620), figure_padding = 8, fontsize = 18)
    cmap = (Reverse(:matter), 1)

    ηm, ηi, ξm, ξi, ri, r2, t, α, sc, ε̇, γ̇, ζ̇, ε̇zz = base
    ncols = 2 * length(ξvals)   # [axis | colorbar] per panel

    grow = 2
    for (col, ξ) in enumerate(ξvals)
        base_col = 2col - 1              # [axis | colorbar]

        params    = merge(base, (ξm = ξ, ξi = ξ))
        _, Pa, _  = analytics_stag(f_anal, X, params, geometry)
        _, Psch,_ = analytics_stag(analytics_schmid_ellipse, X, params, geometry)
        maxdiff   = maximum(abs.(Pa .- Psch))

        # relative pad so a (near-)uniform field doesn't collapse the colorrange
        lo, hi = extrema(Pa)
        pad    = max(abs(hi), 1.0) * 0.05
        crange = (hi - lo) < pad ? (lo - pad, hi + pad) : (lo, hi)

        ax = Axis(fig[grow, base_col], aspect = DataAspect(),
                  title = coltitle(ξ), titlesize = 30,
                  xlabel = L"x", ylabel = L"y", xlabelsize = 28, ylabelsize = 28,
                  xticklabelsize = 18, yticklabelsize = 18)
        hm = heatmap!(ax, xc_p, yc_p, Pa; colorrange = crange, colormap = cmap)
        # only the last (right-most) column carries the p label ("" avoids the
        # empty-LaTeX-string render crash in Makie)
        plabel = (col == length(ξvals)) ? L"p" : ""
        Colorbar(fig[grow, base_col+1], hm, width = 12, label = plabel, labelsize = 34, ticklabelsize = 18)

        text!(ax, 0.5, 0.02;
              text = L"\max|\Delta p|_{\mathrm{incomp}}=%$(fmtnum(maxdiff))",
              space = :relative, align = (:center, :bottom),
              color = :black, fontsize = 22, font = :bold)

        col > 1 && hideydecorations!(ax, grid = false)
    end

    # fixed-parameter caption + clarifying subtitle (ξ is swept, so not listed here);
    # placed after the loop so fig[1,:]/[2,:] actually span every panel column
    Label(fig[1, 1:ncols],
          L"\mu_h=%$(fmtnum(ηm))\;\;\mu_i=%$(fmtnum(ηi))\;\;\varepsilon=%$(fmtnum(ε̇))\;\;\gamma=%$(fmtnum(γ̇))\;\;\zeta=%$(fmtnum(ζ̇))\;\;\alpha=%$(fmtnum(α))\;\;t=%$(fmtnum(t))",
          fontsize = 25, color = :gray25, halign = :center)

    # equalize the four heatmap columns so retained decorations (y-ticks/ylabel on
    # col 1) can't make some panels larger than others.
    for c in (1, 3, 5, 7)
        colsize!(fig.layout, c, Relative(0.2))
    end

    # columns: (ax|cbar) x4 -> keep each colorbar hugging its axis,
    # and put the breathing room between panel groups instead
    colgap!(fig.layout, 8)
    for g in (1, 3, 5, 7);  colgap!(fig.layout, g, 4);  end    # axis <-> its colorbar (tight)
    for g in (2, 4, 6);     colgap!(fig.layout, g, 22); end    # gap between panel groups

    rowgap!(fig.layout, 8)
    rowgap!(fig.layout, 1, 4)   # caption sits close to the panels

    #display(fig)

    mkpath("figures")
    fname = "./figures/Effect_of_compressibility.png"
    save(fname, fig, px_per_unit = 2)
    println("Saved: $fname")
    fig
end
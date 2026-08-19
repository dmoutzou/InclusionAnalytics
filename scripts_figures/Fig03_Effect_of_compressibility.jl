# Include packages and scripts
using ExactFieldSolutions
include("../src/InclusionAnalytics.jl")


define_params(params) = (
    ηm=params.ηm, ηi=params.ηi, ξm=params.ξm, ξi=params.ξi,
    ri=params.r2, t=params.t, α=params.α,
    ε̇=params.ε̇, γ̇=params.γ̇, ζ̇=params.ζ̇, ε̇zz=params.ε̇zz,
)

define_params_schmid_circle(params) = (
    mm=params.ηm, mc=params.ηi, rc=params.r1, gr=params.γ̇, er=params.ε̇,
)

analytics_circle(X; params)  = Stokes2D_Moutzouris_circle(X;  params=define_params(params))
analytics_ellipse(X; params) = Stokes2D_Moutzouris_ellipse(X; params=define_params(params))

set_theme!(theme_latexfonts())

# Schmid & Podladchikov (2003) incompressible elliptical-inclusion solution,
# evaluated with the same interface as analytics_ellipse for direct comparison
function analytics_schmid_circle(X; params)
    x   = @SVector [X[1], X[2]]
    sol = Stokes2D_Schmid2003_circle(x; params=define_params_schmid_circle(params))
    return (V=sol.V, p=sol.p, τ=sol.τ, τzz=0.0)
end

function analytics_schmid_ellipse(X; params)
    x   = @SVector [X[1], X[2]]
    sol = Stokes2D_Schmid2003_ellipse(x; params=define_params(params))
    return (V=sol.V, p=sol.p, τ=sol.τ, τzz=0.0)
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
    f_schmid = geometry == :circular ? analytics_schmid_circle : analytics_schmid_ellipse
    nc       = 501

    Lx, Ly = 1.0, 1.0
    X      = make_grid(nc; Lx=Lx, Ly=Ly)

    scale = 2 / Lx
    xc_p, yc_p = scale .* X.xc, scale .* X.yc

    base   = preset(:Moutzouris2026_1_ps, geometry)
    ξvals  = [0.1, 1.0e0, 1.0e3, 1.0e8]
    incomp = 1.0e8                          # threshold to label as ->infty

    coltitle(ξ) = ξ ≥ incomp ? L"\xi_h=\xi_i\to\infty" : L"\xi_h=\xi_i=%$(fmtnum(ξ))"

    fig  = Figure(size = (1850, 620), figure_padding = 8, fontsize = 18)
    cmap = (Reverse(:matter), 1)

    grow = 2

    Pfields  = Vector{Any}(undef, length(ξvals))
    maxdiffs = Vector{Float64}(undef, length(ξvals))
    for (col, ξ) in enumerate(ξvals)
        params     = merge(base, (ξm = ξ, ξi = ξ))
        _, Pa,  _  = analytics_stag(f_anal,   X, params, geometry)
        _, Psch, _ = analytics_stag(f_schmid, X, params, geometry)
        Pfields[col]  = Pa
        maxdiffs[col] = maximum(abs.(Pa .- Psch))
    end

    qclip  = 0.01
    gain   = 1.5
    allP   = vcat(vec.(Pfields)...)
    lo, hi = quantile(allP, qclip), quantile(allP, 1 - qclip)
    m      = gain * max(abs(lo), abs(hi))
    crange = (-m, m)

    ncols = length(ξvals) + 1   # 4 axes + 1 shared colorbar

    local hm
    for (col, ξ) in enumerate(ξvals)
        ax = Axis(fig[grow, col], aspect = DataAspect(),
                  title = coltitle(ξ), titlesize = 30,
                  xlabel = L"x", ylabel = L"y", xlabelsize = 28, ylabelsize = 28,
                  xticklabelsize = 18, yticklabelsize = 18)
        hm = heatmap!(ax, xc_p, yc_p, Pfields[col]; colorrange = crange, colormap = cmap)

        text!(ax, 0.5, 0.02;
              text = L"\max|\Delta p|_{\mathrm{incomp}}=%$(fmtnum(maxdiffs[col]))",
              space = :relative, align = (:center, :bottom),
              color = :black, fontsize = 22, font = :bold)

        col > 1 && hideydecorations!(ax, grid = false)
    end

    Colorbar(fig[grow, ncols], hm, width = 12, label = L"p",
             labelsize = 34, ticklabelsize = 18)
    Label(fig[1, 1:ncols],
          L"\mu_h=%$(fmtnum(base.ηm))\;\;\mu_i=%$(fmtnum(base.ηi))\;\;\varepsilon=%$(fmtnum(base.ε̇))\;\;\gamma=%$(fmtnum(base.γ̇))\;\;\zeta=%$(fmtnum(base.ζ̇))\;\;\alpha=%$(fmtnum(base.α))\;\;t=%$(fmtnum(base.t))",
          fontsize = 25, color = :gray25, halign = :center)

    for c in 1:length(ξvals)
        colsize!(fig.layout, c, Relative(0.23))
    end

    colgap!(fig.layout, 30)                    # between panels
    colgap!(fig.layout, ncols - 1, 18)

    rowgap!(fig.layout, 8)
    rowgap!(fig.layout, 1, 4)

    mkpath("figures")
    fname = "./figures/Fig03_Effect_of_compressibility.png"
    save(fname, fig, px_per_unit = 2)
    println("Saved: $fname")
    fig
end
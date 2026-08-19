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

analytics_circle(X;  params) = Stokes2D_Moutzouris_circle(X;  params=define_params(params))
analytics_ellipse(X; params) = Stokes2D_Moutzouris_ellipse(X; params=define_params(params))

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
    nc = 501

    Lx, Ly = 1.0, 1.0
    scale  = 2 / Lx

    #circle
    geom_c   = :circular
    params_c = preset(:Schmid2003, geom_c)
    Xc       = make_grid(nc; Lx=Lx, Ly=Ly)
    xc_c, yc_c = scale .* Xc.xc, scale .* Xc.yc

    _, Pnew_circ, _ = analytics_stag(analytics_circle,        Xc, params_c, geom_c)
    _, Psch_circ, _ = analytics_stag(analytics_schmid_circle, Xc, params_c, geom_c)

    #ellipse
    geom_e   = :elliptical
    params_e = preset(:Schmid2003, geom_e)
    Xe       = make_grid(nc; Lx=Lx, Ly=Ly)
    xc_e, yc_e = scale .* Xe.xc, scale .* Xe.yc

    _, Pnew_ell, _ = analytics_stag(analytics_ellipse,        Xe, params_e, geom_e)
    _, Psch_ell, _ = analytics_stag(analytics_schmid_ellipse, Xe, params_e, geom_e)

    #Figure
    set_theme!(theme_latexfonts())
    cmap = Reverse(:matter)
    fig  = Figure(size = (1700, 950), fontsize = 24)

    # layout:  col1 Schmid | col2 This study | col3 shared p colorbar | col4 diff | col5 diff colorbar
    hdrsize = 34
    Label(fig[0, 1], L"\mathrm{Schmid\ \&\ Podladchikov,\ 2003}", fontsize = hdrsize, tellwidth = false)
    Label(fig[0, 2], L"\mathrm{This\ study}",                     fontsize = hdrsize, tellwidth = false)
    Label(fig[0, 4], L"\mathrm{Difference}",                      fontsize = hdrsize, tellwidth = false)

    rows = [
        (Psch_circ, Pnew_circ, xc_c, yc_c),
        (Psch_ell,  Pnew_ell,  xc_e, yc_e),
    ]

    for (row, (ref, new, xx, yy)) in enumerate(rows)
        clim = extrema(ref)
        diff = log10.(abs.(ref .- new))

        bottom = (row == length(rows))
        ax1 = Axis(fig[row,1], aspect=DataAspect(), ylabel=L"y", ylabelsize=34,
                   xlabel=L"x", xlabelsize=34, xticklabelsize=20, yticklabelsize=20)
        ax2 = Axis(fig[row,2], aspect=DataAspect(), xlabel=L"x", xlabelsize=34,
                   xticklabelsize=20, yticklabelsize=20)
        ax3 = Axis(fig[row,4], aspect=DataAspect(), xlabel=L"x", xlabelsize=34,
                   xticklabelsize=20, yticklabelsize=20)

        hm1 = heatmap!(ax1, xx, yy, ref;  colorrange=clim, colormap=cmap)
        hm2 = heatmap!(ax2, xx, yy, new;  colorrange=clim, colormap=cmap)
        hm3 = heatmap!(ax3, xx, yy, diff;                  colormap=cmap)

        # shared colorbar
        Colorbar(fig[row,3], hm1, label = L"p", labelsize = 44, ticklabelsize = 20)
        Colorbar(fig[row,5], hm3, label = L"\log_{10}\,|\Delta p|", labelsize = 31, ticklabelsize = 20)

        # parameters inside the first plot
        if row == 1
            plist = [
                L"\mathrm{\mu_h = %$(fmtnum(params_e.ηm))}",
                L"\mathrm{\mu_i = %$(fmtnum(params_e.ηi))}",
                L"\mathrm{\xi_h = %$(fmtnum(params_e.ξm))}",
                L"\mathrm{\xi_i = %$(fmtnum(params_e.ξi))}",
                L"\mathrm{\varepsilon = %$(fmtnum(params_e.ε̇))}",
                L"\mathrm{\gamma = %$(fmtnum(params_e.γ̇))}",
                L"\mathrm{\zeta = %$(fmtnum(params_e.ζ̇))}",
                L"\mathrm{\alpha = %$(fmtnum(params_e.α))}",
                L"\mathrm{t = %$(fmtnum(params_e.t))}",
            ]
            col1, col2 = plist[1:5], plist[6:end]
            dy = 0.075
            for (k, s) in enumerate(col1)
                text!(ax1, 0.03, 0.03 + (length(col1)-k)*dy; text=s, space=:relative,
                      align=(:left,:bottom), color=:white, fontsize=22)
            end
            for (k, s) in enumerate(col2)
                text!(ax1, 0.34, 0.03 + (length(col2)-k)*dy; text=s, space=:relative,
                      align=(:left,:bottom), color=:white, fontsize=22)
            end
        end

        hideydecorations!(ax2, grid=false)
        hideydecorations!(ax3, grid=false)
        if !bottom
            hidexdecorations!.((ax1, ax2, ax3), grid=false)
        end

        @printf("%-11s  max|Δp| = %.3e   (p range = [%.3e, %.3e])\n",
                row == 1 ? "Circular" : "Elliptical", maximum(abs.(ref .- new)), clim...)
    end

    colgap!(fig.layout, 10)
    rowgap!(fig.layout, 18)
    colgap!(fig.layout, 3, 28)

    save("figures/Fig02_SchmidVsMoutzouris.png", fig; px_per_unit = 200/96)
    fig
end
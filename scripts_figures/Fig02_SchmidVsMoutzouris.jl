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

function analytics_schmid_ellipse(X; params)
    ηm, ηi, ξm, ξi, ri, r2, t, α, sc, ε̇, γ̇, ζ̇, ε̇zz = params
    mc  = ηi / ηm
    Ζ   = to_zeta(sc .* X)
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
    nc = 501

    Lx, Ly = 1.0, 1.0
    scale  = 2 / Lx

    #circle
    geom_c   = :circular
    params_c = preset(:Schmid2003, geom_c)
    Xc       = make_grid(nc; Lx=Lx, Ly=Ly)
    xc_c, yc_c = scale .* Xc.xc, scale .* Xc.yc

    _, Pnew_circ, _ = analytics_stag(analytics_circle, Xc, params_c, geom_c)

    #Schmid (ExactFieldSolutions)
    ηm, ηi, ξm, ξi, ri, r2, t, α, sc, ε̇, γ̇, ζ̇, ε̇zz = params_c
    ps_c = (mm=ηm, mc=ηi, rc=ri, gr=γ̇, er=ε̇)
    Psch_circ = similar(Pnew_circ)
    for j in axes(Psch_circ,2), i in axes(Psch_circ,1)
        Psch_circ[i,j] = Stokes2D_Schmid2003(@SVector([Xc.xc[i], Xc.yc[j]]); params=ps_c).p
    end

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
            pηm, pηi, pξm, pξi, pri, pr2, pt, pα, psc, pε̇, pγ̇, pζ̇, pε̇zz = params_e
            plist = [
                L"\mathrm{\mu_h = %$(fmtnum(pηm))}",
                L"\mathrm{\mu_i = %$(fmtnum(pηi))}",
                L"\mathrm{\xi_h = %$(fmtnum(pξm))}",
                L"\mathrm{\xi_i = %$(fmtnum(pξi))}",
                L"\mathrm{\varepsilon = %$(fmtnum(pε̇))}",
                L"\mathrm{\gamma = %$(fmtnum(pγ̇))}",
                L"\mathrm{\zeta = %$(fmtnum(pζ̇))}",
                L"\mathrm{\alpha = %$(fmtnum(pα))}",
                L"\mathrm{t = %$(fmtnum(pt))}",
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
    rowgap!(fig.layout, 10)
    colgap!(fig.layout, 3, 28)   

    save("figures/Fig02_SchmidVsMoutzouris.png", fig; px_per_unit = 200/96)
    fig
end
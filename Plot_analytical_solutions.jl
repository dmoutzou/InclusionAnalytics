using InclusionAnalytics
using CairoMakie, MathTeXEngine, StaticArrays, Statistics

set_theme!(theme_latexfonts())

function make_grid(nc; Lx=1.0, Ly=1.0)
    Δx, Δy = Lx/nc, Ly/nc
    xce = LinRange(-Lx/2-Δx/2, Lx/2+Δx/2, nc+2); yce = LinRange(-Ly/2-Δy/2, Ly/2+Δy/2, nc+2)
    xc  = LinRange(-Lx/2+Δx/2, Lx/2-Δx/2, nc);   yc  = LinRange(-Ly/2+Δy/2, Ly/2-Δy/2, nc)
    xv  = LinRange(-Lx/2, Lx/2, nc+1);           yv  = LinRange(-Ly/2, Ly/2, nc+1)
    return (xv=xv, yv=yv, xce=xce, yce=yce, xc=xc, yc=yc)
end

# co-located velocity samples on a coarse grid, for quiver overlays
function quiver_arrays(f, params; nq=15, L=0.46)
    xs = LinRange(-L, L, nq); ys = LinRange(-L, L, nq)
    px = Float64[]; py = Float64[]; pu = Float64[]; pv = Float64[]
    for y in ys, x in xs
        vel = f(@SVector([x, y]); params=params).V
        push!(px, x); push!(py, y); push!(pu, vel[1]); push!(pv, vel[2])
    end
    return px, py, pu, pv
end

# two compact lines summarising the parameter set
function param_lines(params)
    ηm, ηi, ξm, ξi, ri, r2, t, α, sc, ε̇, γ̇, ζ̇, ε̇zz = params
    return (
        L"\mu_h=%$(fmtnum(ηm))\;\;\mu_i=%$(fmtnum(ηi))\;\;\xi_h=%$(fmtnum(ξm))\;\;\xi_i=%$(fmtnum(ξi))",
        L"\varepsilon=%$(fmtnum(ε̇))\;\;\gamma=%$(fmtnum(γ̇))\;\;\zeta=%$(fmtnum(ζ̇))\;\;\alpha=%$(fmtnum(α))\;\;t=%$(fmtnum(t))",
    )
end

let
    geometry = :circular
    f_anal   = analytics_circle
    nc       = 501

    Lx, Ly = 1.0, 1.0
    X      = make_grid(nc; Lx=Lx, Ly=Ly)

    # cosmetic axis rescale to [-1, 1] (physics untouched)
    scale = 2 / Lx
    xc_p, yc_p = scale .* X.xc, scale .* X.yc

    # four Duretz forcing cases
    cases = [
        (:Duretz2026_1_ps,   L"\mathrm{Pure\ shear}"),
        (:Duretz2026_2_ss,   L"\mathrm{Simple\ shear}"),
        (:Duretz2026_3_exp,  L"\mathrm{Expansion}"),
        (:Duretz2026_4_comp, L"\mathrm{Compaction}"),
    ]

    fig  = Figure(size = (1300, 1250), fontsize = 18)
    cmap = (Reverse(:matter), 1)

    # grid rows: (header, plot) per panel row -> rows 1,2 and 3,4
    for (idx, (test, title)) in enumerate(cases)
        r    = (idx - 1) ÷ 2 + 1        # 1,1,2,2
        c    = (idx - 1) %  2 + 1        # 1,2,1,2
        base = (c - 1) * 2              # per panel: [axis | colorbar]
        hrow = 2r - 1                    # header row
        prow = 2r                        # plot row

        params   = preset(test, geometry)
        _, Pa, _ = analytics_stag(f_anal, X, params, geometry)
        pl       = param_lines(params)

        # header: title + parameter lines (outside the plot area)
        hg = GridLayout(fig[hrow, base+1])
        Label(hg[1, 1], title, fontsize = 30, tellwidth = false)
        Label(hg[2, 1], pl[1], fontsize = 18, color = :gray25, tellwidth = false)
        Label(hg[3, 1], pl[2], fontsize = 18, color = :gray25, tellwidth = false)
        rowgap!(hg, 3)

        ax = Axis(fig[prow, base+1], aspect = DataAspect(),
                  xlabel = L"x", ylabel = L"y", xlabelsize = 28, ylabelsize = 28,
                  xticklabelsize = 18, yticklabelsize = 18)
        hm = heatmap!(ax, xc_p, yc_p, Pa; colormap = cmap)
        Colorbar(fig[prow, base+2], hm, width = 12, label = L"p", labelsize = 34, ticklabelsize = 18)

        # velocity quivers (length normalised per panel for visibility)
        px, py, pu, pv = quiver_arrays(f_anal, params)
        vmax = maximum(hypot.(pu, pv))
        ls   = vmax == 0 ? 0.0 : 0.8 * (2*0.46*scale/14) / vmax
        acol = test == :Duretz2026_4_comp ? :white : :black
        arrows2d!(ax, scale .* px, scale .* py, pu, pv; lengthscale = ls,
                  tiplength = 7, tipwidth = 7, shaftwidth = 1.0, color = acol)
    end

    colgap!(fig.layout, 8)
    rowgap!(fig.layout, 5)
    rowgap!(fig.layout, 2, 28)   # extra space between the two panel rows
    colgap!(fig.layout, 2, 32)   # gap between the two panel columns

    #display(fig)

    mkpath("figures")
    fname = "./figures/Circular examples.png"
    save(fname, fig, px_per_unit = 2)
    println("Saved: $fname")
end
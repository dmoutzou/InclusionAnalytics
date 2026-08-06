using InclusionAnalytics
using CairoMakie, MathTeXEngine, StaticArrays, Statistics, Printf

# Difference (numerics - analytics) for three forcing cases, side by side.
#   rows    : V_x, V_y, P
#   columns : pure shear, simple shear, mixed
let
    geometry = :elliptical            # :circular | :elliptical
    nc       = 201                    # resolution
    f_anal   = geometry == :circular ? analytics_circle : analytics_ellipse

    cases = [
        (:Duretz2026_1_ps,    L"\mathrm{Pure\ shear}"),
        (:Duretz2026_2_ss,    L"\mathrm{Simple\ shear}"),
        (:Duretz2026_3_exp, L"\mathrm{Expansion}"),
        (:Duretz2026_5_mixed, L"\mathrm{Mixed}"),
    ]

    # --- run every case; keep the error fields, the normalised grid and the params
    results = map(cases) do (test, _)
        params     = preset(test, geometry)
        V, P, σ, X = Stokes2D(nc, test, params, f_anal)
        Va, Pa, σa = analytics_stag(f_anal, X, params, geometry)
        Ve, Pe, σe = errors_stag(V, P, σ, Va, Pa, σa)

        # coordinate normalisation: map the box onto [-1, 1] (physics untouched)
        scale = 2 / (X.xv[end] - X.xv[1])
        g = (xv  = scale .* X.xv,  yv  = scale .* X.yv,
             xc  = scale .* X.xc,  yc  = scale .* X.yc,
             xce = scale .* X.xce, yce = scale .* X.yce)

        return (g = g, fields = (Ve.x, Ve.y, Pe), params = params)
    end

    row_names = [L"u_x", L"u_y", L"P"]

    # staggered grid for each row
    coords(g, row) = row == 1 ? (g.xv,  g.yce) :
                     row == 2 ? (g.xce, g.yv)  :
                                (g.xc,  g.yc)

    # parameter set, split over several lines (one text! call per line)
    param_lines(p) = [
        L"\mu_h=%$(fmtnum(p.ηm))\;\;\mu_i=%$(fmtnum(p.ηi))",
        L"\xi_h=%$(fmtnum(p.ξm))\;\;\xi_i=%$(fmtnum(p.ξi))",
        L"\varepsilon=%$(fmtnum(p.ε̇))\;\;\gamma=%$(fmtnum(p.γ̇))\;\;\zeta=%$(fmtnum(p.ζ̇))",
        L"t=%$(fmtnum(p.t))\;\;\alpha=%$(fmtnum(p.α))",
    ]

    fs     = 23   # base font size
    fs_lbl = 32   # column titles / axis labels
    fs_ann = 14   # parameter annotation

    set_theme!(theme_latexfonts())
    fig = Figure(size = (1500, 1200), fontsize = fs)

    # column titles, spanning [axis | colorbar]
    for (col, (_, title)) in enumerate(cases)
        Label(fig[0, 2col-1:2col], title, tellwidth = false, fontsize = fs_lbl)
    end

    for row in 1:length(row_names), col in 1:length(cases)
        g      = results[col].g
        fld    = results[col].fields[row]
        gx, gy = coords(g, row)

        # ticks exclude the endpoints ±1: those sit on the shared edge between
        # stacked axes and would collide with the neighbouring panel's labels
        ax = Axis(fig[row, 2col-1], aspect = DataAspect(),
                  xlabel = L"x", ylabel = L"y",
                  xlabelsize = fs_lbl, ylabelsize = fs_lbl,
                  xticks = -0.5:0.5:0.5, yticks = -0.5:0.5:0.5,
                  xticklabelsize = 20, yticklabelsize = 20)

        hm = heatmap!(ax, gx, gy, fld; colormap = (Reverse(:matter), 1))

        # parameter set + resolution in the bottom corners of the pressure row only
        if row == length(row_names)
            txt = param_lines(results[col].params)
            x0, y0, dy = -0.95, -0.95, 0.11   # box is [-1, 1] after normalisation
            for (k, s) in enumerate(reverse(txt))
                text!(ax, x0, y0 + (k-1)*dy; text = s, color = :white,
                      align = (:left, :bottom), fontsize = fs_ann)
            end
            text!(ax, 0.95, y0; text = L"n_c = %$(nc)\times%$(nc)", color = :white,
                  align = (:right, :bottom), fontsize = fs_ann)
        end

        # field label only on the rightmost colorbar of each row
        if col == length(cases)
            Colorbar(fig[row, 2col], hm, label = row_names[row],
                     labelsize = fs_lbl, ticklabelsize = 20, width = 14)
        else
            Colorbar(fig[row, 2col], hm, width = 14, ticklabelsize = 20)
        end

        col == 1                || hideydecorations!(ax, grid = false)
        row == length(row_names) || hidexdecorations!(ax, grid = false)
    end

    colgap!(fig.layout, 10)
    rowgap!(fig.layout, 10)

    # Save figure
    mkpath("figures")
    fname = "./figures/differences_$(geometry)_nc$(nc).png"
    save(fname, fig, px_per_unit = 2)
    println("Saved: $fname")
    fig
end
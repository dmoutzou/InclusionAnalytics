# Include packages and scripts
using ExactFieldSolutions
include("../src/InclusionAnalytics.jl")

# Change "geometry" to elliptical to obtain Fig09

define_params(params) = (
    ηm=params.ηm, ηi=params.ηi, ξm=params.ξm, ξi=params.ξi,
    ri=params.r2, t=params.t, α=params.α,
    ε̇=params.ε̇, γ̇=params.γ̇, ζ̇=params.ζ̇, ε̇zz=params.ε̇zz,
)
analytics_circle(X; params)  = Stokes2D_Moutzouris_circle(X;  params=define_params(params))
analytics_ellipse(X; params) = Stokes2D_Moutzouris_ellipse(X; params=define_params(params))

# Convergence study: mean error vs resolution for all four tests, for either
# a circular or an elliptical inclusion.
#Change "geometry" to elliptical to obtain Fig09
let
    geometry = :circular     # :circular | :elliptical
    tests    = [:Moutzouris2026_1_ps, :Moutzouris2026_2_ss, :Moutzouris2026_3_exp, :Moutzouris2026_4_comp]
    titles   = ["Pure Shear", "Simple Shear", "Expansion", "Compaction"]
    nc_list  = [101, 201, 401, 501]
    f_anal   = geometry==:circular ? analytics_circle : analytics_ellipse

    panel = 500  

    okabe_ito = Dict(
        :orange  => RGBf(230/255, 159/255, 0/255),
        :skyblue => RGBf(86/255, 180/255, 233/255),
        :green   => RGBf(0/255, 158/255, 115/255),
        :yellow  => RGBf(240/255, 228/255, 66/255),
        :blue    => RGBf(0/255, 114/255, 178/255),
        :vermil  => RGBf(213/255, 94/255, 0/255),
        :purple  => RGBf(204/255, 121/255, 167/255),
    )

    fields = [
        (:vx,  L"v_x",         okabe_ito[:blue],    :circle),
        (:vy,  L"v_y",         okabe_ito[:orange],  :utriangle),
        (:p,   L"p",           okabe_ito[:green],   :rect),
        (:σxx, L"\sigma_{xx}", okabe_ito[:skyblue], :diamond),
        (:σyy, L"\sigma_{yy}", okabe_ito[:purple],  :dtriangle),
        (:σzz, L"\sigma_{zz}", okabe_ito[:vermil],  :cross),
        (:σxy, L"\sigma_{xy}", okabe_ito[:yellow],  :xcross),
    ]

    # --- LaTeX theme -------------------------------------------------------
    latex_theme = merge(
        theme_latexfonts(),
        Theme(
            fontsize = 22,
            figure_padding = 16,
            Axis = (
                width           = panel,   
                height          = panel,   # the exact size it needs
                titlesize       = 34,
                xlabelsize      = 38,
                ylabelsize      = 38,
                xticklabelsize  = 30,
                yticklabelsize  = 30,
                xminorticksvisible = true,
                yminorticksvisible = true,
                xminorticks     = IntervalsBetween(9),
                yminorticks     = IntervalsBetween(9),
            ),
            Legend = (labelsize = 34, patchsize = (36, 36), framevisible = true,
                      colgap = 18, padding = (12, 12, 8, 8)),
        ),
    )

    fig = with_theme(latex_theme) do

        fig = Figure()
        plots_for_legend  = []
        labels_for_legend = []

        for (k, test) in enumerate(tests)
            params   = preset(test, geometry)
            row, col = divrem(k-1, 2) .+ (1, 1)
            ax = Axis(fig[row, col];
                      xscale = log10, yscale = log10,
                      title  = titles[k],
                      xlabel = L"1/h \;\; ",
                      ylabel = L"\varepsilon")

            invh = Float64[]
            ε   = Dict(sym => Float64[] for (sym, _, _, _) in fields)

            for nc in nc_list

                @info "Test $(test) - nc = $(nc)"

                # Numerics
                V, P, σ, X = Stokes2D(nc, test, params, f_anal)

                # Analytics
                Va, Pa, σa = analytics_stag(f_anal, X, params, geometry)

                # Errors
                Ve, Pe, σe = errors_stag(V, P, σ, Va, Pa, σa)

                Lx = X.xv[end] - X.xv[1]
                push!(invh, nc/Lx)

                push!(ε[:vx],  mean(abs.(Ve.x)))
                push!(ε[:vy],  mean(abs.(Ve.y)))
                push!(ε[:p],   mean(abs.(Pe)))
                push!(ε[:σxx], mean(abs.(σe.xx)))
                push!(ε[:σyy], mean(abs.(σe.yy)))
                push!(ε[:σzz], mean(abs.(σe.zz)))
                push!(ε[:σxy], mean(abs.(σe.xy)))
            end

            for (sym, label, color, marker) in fields
                sc = scatter!(ax, invh, ε[sym];
                    color = color, marker = marker,
                    markersize = 22, strokewidth = 1.4, strokecolor = :black)
                if k == 1
                    push!(plots_for_legend, sc)
                    push!(labels_for_legend, label)
                end
            end

            # O(1) reference line
            xs = [minimum(invh), maximum(invh)]
            y0 = maximum(ε[s][1] for (s, _, _, _) in fields)
            C  = xs[1]*y0*1.3
            ln = lines!(ax, xs, C ./ xs; color = :black, linestyle = :dash, linewidth = 2.5)
            if k == 1
                push!(plots_for_legend, ln)
                push!(labels_for_legend, L"\mathcal{O}(1)")
            end
        end

        Legend(fig[3, 1:2], plots_for_legend, labels_for_legend;
               orientation = :horizontal, nbanks = 1, tellwidth = false)

        rowgap!(fig.layout, 10)
        colgap!(fig.layout, 20)

        resize_to_layout!(fig)

        fig
    end
    if geometry==:circular
        fname = "./figures/Fig08_convergence_$(geometry).png"
        save(fname, fig, px_per_unit = 2)
        display(fig)
        println("Saved: $fname")
    elseif geometry==:elliptical
        fname = "./figures/Fig09_convergence_$(geometry).png"
        save(fname, fig, px_per_unit = 2)
        display(fig)
        println("Saved: $fname")
    end
end
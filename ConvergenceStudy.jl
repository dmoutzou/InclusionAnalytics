using InclusionAnalytics
using CairoMakie, MathTeXEngine, StaticArrays, Statistics

# Convergence study: mean error vs resolution for all four tests, for either
# a circular or an elliptical inclusion.
let
    geometry = :elliptical     # :circular | :elliptical
    tests    = [:Duretz2026_1_ps, :Duretz2026_2_ss, :Duretz2026_3_exp, :Duretz2026_4_comp]
    titles   = ["Pure Shear", "Simple Shear", "Expansion", "Compaction"]
    nc_list  = [101, 201, 401, 501]
    f_anal   = geometry==:circular ? analytics_circle : analytics_ellipse

    panel = 500   # side length in px of each (square) axis

    fields = [
        (:vx,  L"v_x",         :steelblue),
        (:vy,  L"v_y",         :orange),
        (:p,   L"p",           :seagreen),
        (:σxx, L"\sigma_{xx}", :green),
        (:σyy, L"\sigma_{yy}", :yellow),
        (:σzz, L"\sigma_{zz}", :red),
        (:σxy, L"\sigma_{xy}", :black),
    ]

    # --- LaTeX theme -------------------------------------------------------
    latex_theme = merge(
        theme_latexfonts(),
        Theme(
            fontsize = 22,
            figure_padding = 16,
            Axis = (
                width           = panel,   # fixed square axes -> layout knows
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
                      # 1/h = inverse grid spacing (cells per unit length),
                      # with h = L_x/n_c. Log-log, so a slope of -1 is
                      # first-order convergence.
                      xlabel = L"1/h \;\; ",
                      ylabel = L"\varepsilon")

            invh = Float64[]
            ε   = Dict(sym => Float64[] for (sym, _, _) in fields)

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

            for (sym, label, color) in fields
                sc = scatter!(ax, invh, ε[sym]; color = color, markersize = 14)
                if k == 1
                    push!(plots_for_legend, sc)
                    push!(labels_for_legend, label)
                end
            end

            # O(1) reference line
            xs = [minimum(invh), maximum(invh)]
            y0 = maximum(ε[s][1] for (s, _, _) in fields)
            C  = xs[1]*y0*1.3
            ln = lines!(ax, xs, C ./ xs; color = :royalblue, linewidth = 2.5)
            if k == 1
                push!(plots_for_legend, ln)
                push!(labels_for_legend, L"\mathcal{O}(1)")
            end
        end

        Legend(fig[3, 1:2], plots_for_legend, labels_for_legend;
               orientation = :horizontal, nbanks = 1, tellwidth = false)

        rowgap!(fig.layout, 10)
        colgap!(fig.layout, 20)

        # Figure size is derived from the fixed axis sizes -> no white space
        resize_to_layout!(fig)

        fig
    end

    fname = "./figures/convergence_$(geometry).png"
    save(fname, fig, px_per_unit = 2)
    display(fig)
    println("Saved: $fname")
end
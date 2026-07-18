using InclusionAnalytics
using CairoMakie, MathTeXEngine, StaticArrays, Statistics

# Convergence study: mean error vs resolution for all four tests, for either
# a circular or an elliptical inclusion. 
let
    geometry = :elliptical     # :circular | :elliptical
    tests    = [:Schmid2003, :Duretz2026_1, :Duretz2026_2, :Duretz2026_3]
    titles   = ["Test 1", "Test 2", "Test 3", "Test 4"]
    nc_list  = [51, 101, 201, 401]
    f_anal   = geometry==:circular ? analytics_circle : analytics_ellipse

    fields = [
        (:vx,  "vx",   :steelblue),
        (:vy,  "vy",   :orange),
        (:p,   "p",    :seagreen),
        (:σxx, "σxx",  :green),
        (:σyy, "σyy",  :yellow),
        (:σzz, "σzz",  :red),
        (:σxy, "σxy",  :black),
    ]

    fig = Figure(size = (850, 650), fontsize = 12)
    plots_for_legend  = []
    labels_for_legend = String[]

    for (k, test) in enumerate(tests)
        params   = preset(test, geometry)
        row, col = divrem(k-1, 2) .+ (1, 1)
        ax = Axis(fig[row, col]; xscale = log10, yscale = log10,
                   title = titles[k], xlabel = "1/h", ylabel = "ε")

        invh = Float64[]
        ε    = Dict(sym => Float64[] for (sym, _, _) in fields)

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
            sc = scatter!(ax, invh, ε[sym]; color = color, markersize = 9)
            if k == 1
                push!(plots_for_legend, sc)
                push!(labels_for_legend, label)
            end
        end

        # O(1) reference line
        xs = [minimum(invh), maximum(invh)]
        y0 = maximum(ε[sym][1] for (sym, _, _) in fields)
        C  = xs[1]*y0*1.3
        ln = lines!(ax, xs, C ./ xs; color = :royalblue, linewidth = 2)
        if k == 1
            push!(plots_for_legend, ln)
            push!(labels_for_legend, "O(1)")
        end
    end

    Legend(fig[3, 1:2], plots_for_legend, labels_for_legend; orientation = :horizontal, nbanks = 1)

    fname = "./figures/convergence_$(geometry).png"
    save(fname, fig, px_per_unit = 2)
    display(fig)
    println("Saved: $fname")
end

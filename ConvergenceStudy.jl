using InclusionAnalytics
using CairoMakie, MathTeXEngine, StaticArrays, Statistics

# Convergence study: mean error vs resolution for all four tests, for either
# a circular or an elliptical inclusion. 
let
    geometry = :elliptical     # :circular | :elliptical
    a_target = 0.4           # only used when geometry == :elliptical, you can change he size of the inclusion

    tests   = [:Schmid2003, :Duretz2026_1, :Duretz2026_2, :Duretz2026_3]
    titles  = ["Test 1", "Test 2", "Test 3", "Test 4"]
    nc_list = [51, 101, 201, 401]

    fields = [
        (:vx, "vx", :steelblue),
        (:vy, "vy", :orange),
        (:p,  "p",  :seagreen),
    ]

    fig = Figure(size = (850, 650), fontsize = 12)
    plots_for_legend  = []
    labels_for_legend = String[]

    for (k, test) in enumerate(tests)
        row, col = divrem(k-1, 2) .+ (1, 1)
        ax = Axis(fig[row, col]; xscale = log10, yscale = log10,
                   title = titles[k], xlabel = "1/h", ylabel = "ε")

        invh = Float64[]
        ε    = Dict(sym => Float64[] for (sym, _, _) in fields)

        for nc in nc_list
            V, Pt, S, X, params, scale = Stokes2D(nc, test; geometry=geometry, a_target=a_target)
            Lx = X.xv[end] - X.xv[1]
            push!(invh, nc/Lx)

            Vx_an = [f_anal(@SVector([X.xv[i],  X.yce[j]]); params=params, geometry=geometry, S=scale).V[1] for i in axes(V.x, 1), j in axes(V.x, 2)]
            Vy_an = [f_anal(@SVector([X.xce[i], X.yv[j] ]); params=params, geometry=geometry, S=scale).V[2] for i in axes(V.y, 1), j in axes(V.y, 2)]
            P_an  = [f_anal(@SVector([X.xc[i],  X.yc[j] ]); params=params, geometry=geometry, S=scale).p    for i in axes(Pt,  1), j in axes(Pt,  2)]

            push!(ε[:vx], mean(abs.(Vx_an .- V.x)))
            push!(ε[:vy], mean(abs.(Vy_an .- V.y)))
            push!(ε[:p],  mean(abs.(P_an  .- Pt )))
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

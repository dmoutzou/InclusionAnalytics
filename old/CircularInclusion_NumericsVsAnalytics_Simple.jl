using InclusionAnalytics
using CairoMakie, MathTeXEngine
using ExactFieldSolutions, JLD2, StaticArrays, Printf, Statistics

# Compare solutions
let
    test        = :Schmid2003   # :Schmid2003 | :Duretz2026_1 | :Duretz2026_2 | :Duretz2026_3
    formulation = :new          # :old | :new
    nc          = 200           # resolution

    Lx, Ly = 1.0, 1.0
    L1_error, V_num, P_num, σ_num, Va, Pa, Sa, X =
        Stokes2D(nc, test; formulation=formulation)

    Δx = Lx/nc;  Δy = Ly/nc

    #Analytical on the staggered grids (using the same solution that set the BCs)
    f_anal = formulation == :new ? Analytics_new : Analytics_old
    bc_params_full = preset_circle(test)
    params_bc = (ηm=bc_params_full.ηm, ηi=bc_params_full.ηi,
                 ξm=bc_params_full.ξm, ξi=bc_params_full.ξm,
                 γ̇=bc_params_full.γ̇, ε̇=bc_params_full.ε̇, ζ̇=bc_params_full.ζ̇,
                 ri=bc_params_full.ri)

    Vx_an = zeros(nc+1, nc+2);  Vy_an = zeros(nc+2, nc+1);  P_an = zeros(nc, nc)
    for I in CartesianIndices(Vx_an)
        i, j = I[1], I[2]
        Vx_an[i,j] = f_anal(@SVector([X.xv[i]; X.yce[j]]); params=params_bc).V[1]
    end
    for I in CartesianIndices(Vy_an)
        i, j = I[1], I[2]
        Vy_an[i,j] = f_anal(@SVector([X.xce[i]; X.yv[j]]); params=params_bc).V[2]
    end
    for I in CartesianIndices(P_an)
        i, j = I[1], I[2]
        P_an[i,j] = f_anal(@SVector([X.xc[i]; X.yc[j]]); params=params_bc).p
    end

    dVx = Vx_an .- V_num.x
    dVy = Vy_an .- V_num.y
    dP  = P_an  .- P_num 

    # @printf("max|dVx| = %1.3e   max|dVy| = %1.3e   max|dP| = %1.3e\n",
    #         maximum(abs, dVx), maximum(abs, dVy), maximum(abs, dP))
    @printf("mean|dVx| = %1.3e   mean|dVy| = %1.3e   mean|dP| = %1.3e\n",
        mean(abs, dVx), mean(abs, dVy), mean(abs, dP))

    field_data = [
        ("Vx", X.xv,  X.yce, Vx_an, V_num.x, dVx),
        ("Vy", X.xce, X.yv,  Vy_an, V_num.y, dVy),
        ("P",  X.xc,  X.yc,  P_an,  P_num,   dP ),
    ]

    fig = Figure(size=(1100, 900), fontsize=12)
    Label(fig[-1, 1:6], "nc=$nc  |  test=$test  |  formulation=$formulation",
          tellwidth=false, fontsize=15, font=:bold)
    Label(fig[0, 1:2], "Analytical",  tellwidth=false, fontsize=14, font=:bold)
    Label(fig[0, 3:4], "Numerical",   tellwidth=false, fontsize=14, font=:bold)
    Label(fig[0, 5:6], "Difference",  tellwidth=false, fontsize=14, font=:bold)

    for (row, (name, gx, gy, an, num, diff)) in enumerate(field_data)
        an_c  = an  .- mean(an)
        num_c = num .- mean(num)
        dlim  = max(maximum(abs, diff), 1e-12)

        ax1 = Axis(fig[row, 1], aspect=DataAspect(), ylabel=name)
        hm1 = heatmap!(ax1, gx, gy, an_c;  colormap=(Reverse(:matter), 1))
        Colorbar(fig[row, 2], hm1, width=12)

        ax2 = Axis(fig[row, 3], aspect=DataAspect())
        hm2 = heatmap!(ax2, gx, gy, num_c; colormap=(Reverse(:matter), 1))
        Colorbar(fig[row, 4], hm2, width=12)

        ax3 = Axis(fig[row, 5], aspect=DataAspect())
        hm3 = heatmap!(ax3, gx, gy, diff;  colormap=(Reverse(:matter), 1))
        Colorbar(fig[row, 6], hm3, width=12)
    end

    fname = "./figures/comparison_nc$(nc)_$(test)_$(formulation).png"
    save(fname, fig, px_per_unit=2)
    display(fig)
    println("Saved: $fname")
end

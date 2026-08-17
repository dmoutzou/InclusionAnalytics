#Include packages and scripts
using ExactFieldSolutions
include("src/InclusionAnalytics.jl")

#Define parameters
define_params(params) = (
    ηm=params.ηm, ηi=params.ηi, ξm=params.ξm, ξi=params.ξi,
    ri=params.r2, t=params.t, α=params.α,
    ε̇=params.ε̇, γ̇=params.γ̇, ζ̇=params.ζ̇, ε̇zz=params.ε̇zz,
)
analytics_circle(X; params)  = Stokes2D_Moutzouris_circle(X;  params=define_params(params))
analytics_ellipse(X; params) = Stokes2D_Moutzouris_ellipse(X; params=define_params(params))

let
    geometries = [:circular, :elliptical]
    titles     = ["Circular inclusion", "Elliptical inclusion"]

    test      = :Duretz2026_7_oop
    test_ezz0 = :Duretz2026_8_zeroop
    nc_list   = [101, 201, 401, 501]

    # Geometry 
    oop_geometry = :elliptical
    nc_field     = 301

    okabe_ito = Dict(
        :blue    => RGBf(0/255, 114/255, 178/255),
        :orange  => RGBf(230/255, 159/255, 0/255),
        :green   => RGBf(0/255, 158/255, 115/255),
        :vermil  => RGBf(213/255, 94/255, 0/255),
    )

    fields = [
        (:vx,  L"v_x",         okabe_ito[:blue],   :circle),
        (:vy,  L"v_y",         okabe_ito[:orange], :utriangle),
        (:p,   L"p",           okabe_ito[:green],  :rect),
        (:σzz, L"\sigma_{zz}", okabe_ito[:vermil], :diamond),
    ]

    pcmap = Reverse(:matter)

    function panel_label!(ax, label::AbstractString)
        poly!(ax, Rect(0.03, 0.03, 0.10, 0.10);
            space = :relative,
            color = :white,
            strokecolor = :black,
            strokewidth = 1.5,
        )
        text!(ax, 0.08, 0.08;
            text = label,
            space = :relative,
            align = (:center, :center),
            fontsize = 26,
            font = :bold,
            color = :black,
        )
    end

    fmtnum(x) = @sprintf("%.3g", x)

    function param_list!(ax, params)
        # ηm, ηi, ξm, ξi, ri, r2, t, α, sc, ε̇, γ̇, ζ̇, ε̇zz
        ηm, ηi, ξm, ξi, ri, r2, t, α, sc, ε̇, γ̇, ζ̇, ε̇zz = params
        plist = [
            L"\mathrm{\mu_h = %$(fmtnum(ηm))}",
            L"\mathrm{\mu_i = %$(fmtnum(ηi))}",
            L"\mathrm{\xi_h = %$(fmtnum(ξm))}",
            L"\mathrm{\xi_i = %$(fmtnum(ξi))}",
            L"\mathrm{\varepsilon = %$(fmtnum(ε̇))}",
            L"\mathrm{\gamma = %$(fmtnum(γ̇))}",
            L"\mathrm{\zeta = %$(fmtnum(ζ̇))}",
            L"\mathrm{\varepsilon_{zz} = %$(fmtnum(ε̇zz))}",
            L"\mathrm{\alpha = %$(fmtnum(α))}",
            L"\mathrm{t = %$(fmtnum(t))}",
        ]
        col1, col2 = plist[1:5], plist[6:end]
        dy = 0.075
        for (k, s) in enumerate(col1)
            text!(ax, 0.66, 0.03 + (length(col1)-k)*dy; text=s, space=:relative,
                  align=(:left,:bottom), color=:white, fontsize=20)
        end
        for (k, s) in enumerate(col2)
            text!(ax, 0.97, 0.03 + (length(col2)-k)*dy; text=s, space=:relative,
                  align=(:right,:bottom), color=:white, fontsize=20)
        end
    end

    latex_theme = merge(
        theme_latexfonts(),
        Theme(
            fontsize = 22,
            figure_padding = 16,
            Axis = (
                width = 520,
                height = 520,
                titlesize = 34,
                xlabelsize = 38,
                ylabelsize = 38,
                xticklabelsize = 30,
                yticklabelsize = 30,
                xminorticksvisible = true,
                yminorticksvisible = true,
                xminorticks = IntervalsBetween(9),
                yminorticks = IntervalsBetween(9),
            ),
            Legend = (
                labelsize = 34,
                patchsize = (36,36),
                framevisible = true,
                colgap = 18,
            ),
        )
    )

    fig = with_theme(latex_theme) do

        fig = Figure()

        f_anal_field = oop_geometry == :circular ? analytics_circle : analytics_ellipse

        params_ezz0 = preset(test_ezz0, oop_geometry)
        V0, P0, σ0, X0 = Stokes2D(nc_field, test_ezz0, params_ezz0, f_anal_field)

        params_full = preset(test, oop_geometry)
        V1, P1, σ1, X1 = Stokes2D(nc_field, test, params_full, f_anal_field)

        clim = extrema((extrema(P0)..., extrema(P1)...))

        ax1 = Axis(fig[1,1],
            title  = L"\text{Plane-strain, } \varepsilon_{zz} = 0",
            xlabel = L"x",
            ylabel = L"y",
            aspect = DataAspect(),
        )
        hm1 = heatmap!(ax1, X0.xc, X0.yc, P0; colorrange = clim, colormap = pcmap)
        panel_label!(ax1, "A")

        ax2 = Axis(fig[1,2],
            title  = L"\text{Plane-strain, } \varepsilon_{zz} = 1.0",
            xlabel = L"x",
            ylabel = L"y",
            aspect = DataAspect(),
        )
        hm2 = heatmap!(ax2, X1.xc, X1.yc, P1; colorrange = clim, colormap = pcmap)
        panel_label!(ax2, "B")
        param_list!(ax2, params_full)

        hideydecorations!(ax2, grid = false)

        Colorbar(fig[1,3], hm1, label = L"p", labelsize = 38, ticklabelsize = 24)

        linkaxes!(ax1, ax2)

        plots = []
        labels = []
        conv_labels = ["C", "D"]

        for (k, geometry) in enumerate(geometries)

            f_anal = geometry == :circular ?
                analytics_circle :
                analytics_ellipse

            params = preset(test, geometry)

            ax = Axis(fig[2,k],
                xscale = log10,
                yscale = log10,
                title = titles[k],
                xlabel = L"1/h",
                ylabel = L"\varepsilon"
            )

            invh = Float64[]
            ε = Dict(sym => Float64[] for (sym,_,_,_) in fields)

            for nc in nc_list

                @info "$geometry : nc = $nc"

                V,P,σ,X = Stokes2D(nc,test,params,f_anal)

                Va,Pa,σa = analytics_stag(f_anal,X,params,geometry)

                Ve,Pe,σe = errors_stag(V,P,σ,Va,Pa,σa)

                Lx = X.xv[end]-X.xv[1]

                push!(invh,nc/Lx)

                push!(ε[:vx],  mean(abs.(Ve.x)))
                push!(ε[:vy],  mean(abs.(Ve.y)))
                push!(ε[:p],   mean(abs.(Pe)))
                push!(ε[:σzz], mean(abs.(σe.zz)))
            end

            for (sym,label,color,marker) in fields

                sc = scatter!(
                    ax,
                    invh,
                    ε[sym];
                    color = color,
                    marker = marker,
                    markersize = 16,
                    strokewidth = 1.2,
                    strokecolor = :black,
                )

                if k==1
                    push!(plots,sc)
                    push!(labels,label)
                end
            end

            xs = [minimum(invh), maximum(invh)]
            y0 = maximum(ε[s][1] for (s,_,_,_) in fields)
            C = xs[1]*y0*1.3

            ln = lines!(
                ax,
                xs,
                C ./ xs;
                color = :black,
                linestyle = :dash,
                linewidth = 2.5,
            )

            if k==1
                push!(plots,ln)
                push!(labels,L"\mathcal{O}(1)")
            end

            panel_label!(ax, conv_labels[k])
        end

        Legend(fig[3,1:2],
            plots,
            labels,
            orientation=:horizontal,
            tellwidth=false)

        colgap!(fig.layout,20)
        rowgap!(fig.layout,10)

        resize_to_layout!(fig)

        fig
    end

    save("./figures/Fig10_OOP.png",fig,px_per_unit=2)

    display(fig)

end
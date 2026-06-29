using CairoMakie, JLD2, MathTeXEngine
Makie.update_theme!( fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))
Makie.inline!(true)

let 

    f = jldopen("TruncationSystematics.jld2")
    L1 = f["errors"]

    L1.h[1,:]

    fig = Figure( fontsize=18 )

    ax = Axis(fig[1,1], title=L"Test $1$", xscale = log10, yscale = log10, xlabel=L"$1 / h$", ylabel = L"$\epsilon$")
    scatter!(ax, 1 ./ L1.h[1,:],  L1.Vx[1,:])
    scatter!(ax, 1 ./ L1.h[1,:],  L1.Vy[1,:])
    scatter!(ax, 1 ./ L1.h[1,:],   L1.P[1,:])
    scatter!(ax, 1 ./ L1.h[1,:], L1.σxx[1,:])
    scatter!(ax, 1 ./ L1.h[1,:], L1.σyy[1,:])
    scatter!(ax, 1 ./ L1.h[1,:], L1.σzz[1,:])
    scatter!(ax, 1 ./ L1.h[1,:], L1.σxy[1,:])
    lines!(ax, 1 ./ L1.h[1,:], L1.h[1,:])

    ax = Axis(fig[1,2], title=L"Test $2$", xscale = log10, yscale = log10, xlabel=L"$1 / h$", ylabel = L"$\epsilon$")
    scatter!(ax, 1 ./ L1.h[1,:],  L1.Vx[2,:])
    scatter!(ax, 1 ./ L1.h[1,:],  L1.Vy[2,:])
    scatter!(ax, 1 ./ L1.h[1,:],   L1.P[2,:])
    scatter!(ax, 1 ./ L1.h[1,:], L1.σxx[2,:])
    scatter!(ax, 1 ./ L1.h[1,:], L1.σyy[2,:])
    scatter!(ax, 1 ./ L1.h[1,:], L1.σzz[2,:])
    scatter!(ax, 1 ./ L1.h[1,:], L1.σxy[2,:])
    lines!(ax, 1 ./ L1.h[1,:], L1.h[1,:])

    ax = Axis(fig[2,1], title=L"Test $3$", xscale = log10, yscale = log10, xlabel=L"$1 / h$", ylabel = L"$\epsilon$")
    scatter!(ax, 1 ./ L1.h[1,:],  L1.Vx[3,:])
    scatter!(ax, 1 ./ L1.h[1,:],  L1.Vy[3,:])
    scatter!(ax, 1 ./ L1.h[1,:],   L1.P[3,:])
    scatter!(ax, 1 ./ L1.h[1,:], L1.σxx[3,:])
    scatter!(ax, 1 ./ L1.h[1,:], L1.σyy[3,:])
    scatter!(ax, 1 ./ L1.h[1,:], L1.σzz[3,:])
    scatter!(ax, 1 ./ L1.h[1,:], L1.σxy[3,:])
    lines!(ax, 1 ./ L1.h[1,:], L1.h[1,:])

    ax = Axis(fig[2,2], title=L"Test $4$", xscale = log10, yscale = log10, xlabel=L"$1 / h$", ylabel = L"$\epsilon$")
    scatter!(ax, 1 ./ L1.h[1,:],  L1.Vx[4,:], label=L"$v_x$")
    scatter!(ax, 1 ./ L1.h[1,:],  L1.Vy[4,:], label=L"$v_y$")
    scatter!(ax, 1 ./ L1.h[1,:],   L1.P[4,:], label=L"$p$")
    scatter!(ax, 1 ./ L1.h[1,:], L1.σxx[4,:], label=L"$\sigma_{xx}$")
    scatter!(ax, 1 ./ L1.h[1,:], L1.σyy[4,:], label=L"$\sigma_{yy}$")
    scatter!(ax, 1 ./ L1.h[1,:], L1.σzz[4,:], label=L"$\sigma_{zz}$")
    scatter!(ax, 1 ./ L1.h[1,:], L1.σxy[4,:], label=L"$\sigma_{xy}$")
    lines!(ax, 1 ./ L1.h[1,:], L1.h[1,:], label=L"$O(1)$")

    Legend(fig[3, 1:2], ax;
    framevisible = false,
    orientation = :horizontal,
    tellwidth = false,
    tellheight = true)

    display(fig)
    
end
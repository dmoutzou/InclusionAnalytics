using CairoMakie, MathTeXEngine, Printf, Statistics, LinearAlgebra, ExactFieldSolutions, ForwardDiff, StaticArrays
using BenchmarkTools

# ---- Matrix potentials ----
ϕ_mat(z,   k, γ̇, ε̇, ζ̇, ηm, κm, ri) =  2*ηm*ζ̇/(κm-1)*z - im/2*ηm*γ̇*z      -   (im*γ̇ + 2*ε̇) * ηm * k[1] * ri^2 / z^1 
ϕ′_mat(z,  k, γ̇, ε̇, ζ̇, ηm, κm, ri) =  2*ηm*ζ̇/(κm-1)   - im/2*ηm*γ̇        +   (im*γ̇ + 2*ε̇) * ηm * k[1] * ri^2 / z^2 
ϕ′′_mat(z, k, γ̇, ε̇, ζ̇, ηm, κm, ri) =                  - 2*(im*γ̇ + 2*ε̇) * ηm * k[1] * ri^2 / z^3       
ψ_mat(z,   k, γ̇, ε̇, ζ̇, ηm, κm, ri) =                   (im*γ̇ - 2*ε̇)*ηm*z -   (im*γ̇ + 2*ε̇) * ηm * k[1] * ri^4 / z^3       
ψ′_mat(z,  k, γ̇, ε̇, ζ̇, ηm, κm, ri) =                   (im*γ̇ - 2*ε̇)*ηm   + 3*(im*γ̇ + 2*ε̇) * ηm * k[1] * ri^4 / z^4       

# ---- Inclusion potentials ----
ϕ_inc(z,   k, γ̇, ε̇, ζ̇, ηi, κi, ri) =  2*ηi*k[4]*ζ̇/(κi-1)*z  - im/2*ηi*γ̇*z       
ϕ′_inc(z,  k, γ̇, ε̇, ζ̇, ηi, κi, ri) =  2*ηi*k[4]*ζ̇/(κi-1)    - im/2*ηi*γ̇  
ϕ′′_inc(z, k, γ̇, ε̇, ζ̇, ηi, κi, ri) =  0.0
ψ_inc(z,   k, γ̇, ε̇, ζ̇, ηi, κi, ri) =  2*(im*γ̇ - 2*ε̇) * ηi * k[2] * z
ψ′_inc(z,  k, γ̇, ε̇, ζ̇, ηi, κi, ri) =  2*(im*γ̇ - 2*ε̇) * ηi * k[2]    

function constants_analytical(ri, ηm, ηi, ξm, ξi, γ̇, ε̇, ζ̇)
    νm = (3*ξm - 2*ηm) / (2*(3*ξm + ηm))
    νi = (3*ξi - 2*ηi) / (2*(3*ξi + ηi))
    κm = 3 - 4*νm
    κi = 3 - 4*νi
    Mm = ξm + ηm/3
    Mi = ξi + ηi/3
    # ---- Deviatoric constants ----
    k1 = (ηi - ηm) / (ηm + ηi*κm)
    k2 = ηm*(κm + 1) / (2*(ηm + ηi*κm))
    # ---- Volumetric constants ----
    k3 = (Mm - Mi) / (Mi + ηm)   # matrix perturbation amplitude / ζ̇
    k4 = 1 + k3                  # inclusion volumetric strain rate / ζ̇
    return @SVector([k1, k2, k3, k4])         
end

function Analytics(X; params=params)
    ηm, ηi, ξm, ξi, ri, γ̇, ε̇, ζ̇ = params
    z   = X[1] + im*X[2]
    z̄   = conj(z)
    r   = abs(z)
    θ   = angle(z)

    νm  = (3*ξm - 2*ηm)/(2*(3*ξm + ηm))
    νi  = (3*ξi - 2*ηi)/(2*(3*ξi + ηi))
    κm  = 3 - 4*νm
    κi  = 3 - 4*νi

    # ---- Analytical n=0 constants ----
    k = constants_analytical(ri, ηm, ηi, ξm, ξi, γ̇, ε̇, ζ̇)

    if r > ri
        η, κ, ν = ηm, κm, νm

        ϕ′  = ϕ′_mat(z,  k, γ̇, ε̇, ζ̇, ηm, κm, ri)
        ϕ′′ = ϕ′′_mat(z, k, γ̇, ε̇, ζ̇, ηm, κm, ri)
        ϕ   = ϕ_mat(z,   k, γ̇, ε̇, ζ̇, ηm, κm, ri)
        ψ   = ψ_mat(z,   k, γ̇, ε̇, ζ̇, ηm, κm, ri)
        ψ′  = ψ′_mat(z,  k, γ̇, ε̇, ζ̇, ηm, κm, ri)

        # n=0 perturbation velocity: vr = k_pert*rc²/r → W = k_pert*rc²/z̄
        W_pert   =  k[3] * ζ̇ * ri^2/z̄

        # n=0 perturbation stress (divergence-free, purely deviatoric)
        # σrr_pert = -2ηm*k_pert*rc²/r²
        # σθθ_pert = +2ηm*k_pert*rc²/r²
        # σrθ_pert = 0
        fac      =  2 * ηm * k[3] * ζ̇ * ri^2/r^2
        sxx_pert = -fac*cos(2θ)              # = (σrr*cos²θ + σθθ*sin²θ)
        syy_pert =  fac*cos(2θ)              # = (σrr*sin²θ + σθθ*cos²θ)
        sxy_pert = -fac*sin(2θ)              # = (σrr-σθθ)*sinθcosθ

    else
        η, κ, ν = ηi, κi, νi
        ϕ′  = ϕ′_inc(z,  k, γ̇, ε̇, ζ̇, ηi, κi, ri)
        ϕ′′ = ϕ′′_inc(z, k, γ̇, ε̇, ζ̇, ηi, κi, ri)
        ϕ   = ϕ_inc(z,   k, γ̇, ε̇, ζ̇, ηi, κi, ri)
        ψ   = ψ_inc(z,   k, γ̇, ε̇, ζ̇, ηi, κi, ri)
        ψ′  = ψ′_inc(z,  k, γ̇, ε̇, ζ̇, ηi, κi, ri)

        W_pert   = 0.0 + 0.0im
        sxx_pert = 0.0
        syy_pert = 0.0
        sxy_pert = 0.0
    end

    # ---- Stress from potentials (far-field / inclusion uniform field) ----
    Q    = z̄*ϕ′′ + ψ′                        # = 0 for both (ϕ''=ψ'=0)
    sxx  = 2*real(ϕ′) - real(Q) + sxx_pert
    syy  = 2*real(ϕ′) + real(Q) + syy_pert
    sxy  = imag(Q)              + sxy_pert
    szz  = ν*(sxx + syy)
    p    = -1/3*(sxx + syy + szz)

    # ---- Velocity from Muskhelishvili + n=0 perturbation ----
    W    = 1/(2*η)*(κ*ϕ - z*conj(ϕ′) - conj(ψ)) + W_pert

    return (V  = @SVector([real(W), imag(W)]),
            p  = p,
            τ  = @SMatrix([sxx+p  sxy; sxy  syy+p]),
            τzz= szz + p)
end

let 
    # Parameters
    ηm = 1.0
    ηi = 1e-3
    ξm = 1e0
    ξi = 1e0
    ri = 0.1
    t  = π/3
    γ̇  = 1.0
    ε̇  = 1e-10
    ζ̇  = 0.0
    params = (ηm=ηm, ηi=ηi, ξm=ξm, ξi=ξi, ri=ri, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇)

    X = @SVector([0.0, 0.0])
    @benchmark sol = Analytics($X; params=$params )

    # Constants
    # k = constants(ri, ηm, ηi, ξm, ξi, γ̇, ε̇, ζ̇)
    # # @show k
    # k = constants_analytical(ri, ηm, ηi, ξm, ξi, γ̇, ε̇, ζ̇)
    # @show k
end


let 
    # test = :Schmid2003 
    # test = :Duretz2026_1   
    test = :Duretz2026_2

    n  = (x=101, y=101)
    x  = LinRange(-1/2, 1/2, n.x)
    y  = LinRange(-1/2, 1/2, n.y)
    ri = 0.1
    if test == :Schmid2003
        ηm, ηi = 1.0, 1e-2
        ξm, ξi = 1e10, 1e10
        ε̇      =-1e0
        γ̇      = 1.0
        ζ̇      = 0.0
    elseif test==:Duretz2026_1
        ηm, ηi = 1.0, 1e-2
        ξm, ξi = 1e0, 1e0
        ε̇      =-1e0
        γ̇      = 0.0
        ζ̇      = 0.0
    elseif test==:Duretz2026_2
        ηm, ηi = 1.0, 1e-1
        ξm, ξi = 1e0, 1e0
        ε̇      = 1e-10
        γ̇      = 0.0
        ζ̇      = 1e0
    end
    params            = (ηm=ηm, ηi=ηi, ξm=ξm, ξi=ξi, ri=ri, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇)
    params_Schmid2003 = (mm=ηm, mc=ηi, rc=ri, gr=γ̇, er=ε̇)

    V = (
        x = zeros(n...),
        y = zeros(n...)
    )
    σ = (
        xx = zeros(n...),
        xy = zeros(n...)
    )
    P = zeros(n...)

    V_Schmid2003 = (
        x = zeros(n...),
        y = zeros(n...)
    )
    σ_Schmid2003 = (
        xx = zeros(n...),
        xy = zeros(n...)
    )
    P_Schmid2003 = zeros(n...)
    
    @time for I in CartesianIndices(V.x)
        i, j = I[1], I[2]
        X = @SVector([x[i]; y[j]])
        # Call test
        sol = Analytics( X; params=params )
        V.x[i,j]  = sol.V[1]
        V.y[i,j]  = sol.V[2]
        P[i,j]    = sol.p 
        σ.xx[i,j] = sol.τ[1,1] - sol.p
        σ.xy[i,j] = sol.τ[1,2] 
        # Call ExactFieldSolutions
        sol = Stokes2D_Schmid2003( X; params=params_Schmid2003 )
        V_Schmid2003.x[i,j]  = sol.V[1]
        V_Schmid2003.y[i,j]  = sol.V[2]
        P_Schmid2003[i,j]    = sol.p  
        σ_Schmid2003.xx[i,j] = sol.τ[1,1] - sol.p 
        σ_Schmid2003.xy[i,j] = sol.τ[1,2]
    end

    # Checks
    diff = V.y .+  V.x'
    @show extrema(diff)
    @show extrema(V_Schmid2003.x .- V.x)
    @show extrema(V_Schmid2003.y .- V.y)
    @show extrema(P_Schmid2003 .- P)

    if test == :Schmid2003
        col_p   = (-3,3)
        col_sxx = (-5,0)
        col_sxy = (-0.75,0.75)
    else    
        col_p   = extrema(P)
        col_sxx = extrema(σ.xx)
        col_sxy = extrema(σ.xy)
    end

    # Visualise
    fig = Figure(size=(900, 1000))

    ax  = Axis(fig[1,1], aspect=DataAspect(), title=L"$$This study")
    hm  = heatmap!(ax, x, y, V.x, colormap=(Reverse(:matter), 1))
    st  = 10
    arrows2d!(ax, x[1:st:end], y[1:st:end], V.x[1:st:end,1:st:end], V.y[1:st:end,1:st:end], lengthscale=0.3 )
    Colorbar(fig[1,2], hm)
    ax  = Axis(fig[2,1], aspect=DataAspect())
    hm  = heatmap!(ax, x, y, V.y, colormap=(Reverse(:matter), 1))
    Colorbar(fig[2,2], hm)
    ax  = Axis(fig[3,1], aspect=DataAspect())
    hm  = heatmap!(ax, x, y, P, colorrange=col_p, colormap=(Reverse(:matter), 1))
    Colorbar(fig[3,2], hm)
    ax  = Axis(fig[4,1], aspect=DataAspect())
    hm  = heatmap!(ax, x, y, σ.xx, colorrange=col_sxx, colormap=(Reverse(:matter), 1))
    Colorbar(fig[4,2], hm)
    ax  = Axis(fig[5,1], aspect=DataAspect())
    hm  = heatmap!(ax, x, y, σ.xy, colorrange=col_sxy, colormap=(Reverse(:matter), 1))
    Colorbar(fig[5,2], hm)

    ax  = Axis(fig[1,3], aspect=DataAspect(), title=L"$$Schmid & Podladchikov (2004)")
    hm  = heatmap!(ax, x, y, V_Schmid2003.x .- 0.0*V.x, colormap=(Reverse(:matter), 1))
    Colorbar(fig[1,4], hm)
    ax  = Axis(fig[2,3], aspect=DataAspect())
    hm  = heatmap!(ax, x, y, V_Schmid2003.y .- 0*V.y, colormap=(Reverse(:matter), 1))
    Colorbar(fig[2,4], hm)
    ax  = Axis(fig[3,3], aspect=DataAspect())
    hm  = heatmap!(ax, x, y, P_Schmid2003 .-0*P, colorrange=(-3,3), colormap=(Reverse(:matter), 1))
    Colorbar(fig[3,4], hm)
    ax  = Axis(fig[4,3], aspect=DataAspect())
    hm  = heatmap!(ax, x, y, σ_Schmid2003.xx .-0*σ.xx, colorrange=(-5,0), colormap=(Reverse(:matter), 1))
    Colorbar(fig[4,4], hm)
    ax  = Axis(fig[5,3], aspect=DataAspect())
    hm  = heatmap!(ax, x, y, σ_Schmid2003.xy .-0*σ.xy, colorrange=(-0.75,0.75), colormap=(Reverse(:matter), 1))
    Colorbar(fig[5,4], hm)
    display(fig)

    @show params

end


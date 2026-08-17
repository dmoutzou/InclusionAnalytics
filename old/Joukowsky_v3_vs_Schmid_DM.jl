using StaticArrays, CairoMakie

# physical position from the conformal coordinate
joukowski(λ) = λ + 1.0/λ          # ω(λ) = z   
t_to_ri(t) = sqrt((t - 1.0)*(t + 1.0)) / (t - 1.0)

function Analytics_new(X; params)
    ηm, ηi, ξm, ξi, t, γ̇, ε̇, ζ̇ , α = params
    ri      = t_to_ri(t)
    τ, σ    = X[1], X[2]
    #Kolosov
    νm      = (3.0ξm - 2.0ηm) / (2.0*(3.0ξm + ηm))
    νi      = (3.0ξi - 2.0ηi) / (2.0*(3.0ξi + ηi))
    κm      = 3 - 4νm
    κi      = 3 - 4νi
    #Boundary terms with angle of rotation
    P_m     = -im*ηm*γ̇/(κm+1.0) + 2.0*ηm*ζ̇/(κm-1)
    Q_m     = (im*γ̇ - 2.0ε̇)*ηm*exp(-2.0*im*α)
    P_mR    = real(P_m); P_mI = imag(P_m)
    Q_mR    = real(Q_m); Q_mI = imag(Q_m)
    #radius factor
    r4      = ri^4
    #Solution coefficients
    K       = -2.0*ηi*κm*P_mR + 2.0*ηi*P_mR + 2.0*ηi*Q_mR*ri^2 + 2.0*ηm*κi*P_mR +
                ηm*κi*Q_mR*ri^2 - 2.0*ηm*P_mR - ηm*Q_mR*ri^2
    L       = ηi*κm*P_mR*r4 - ηi*P_mR - ηi*Q_mR*ri^2 + ηm*P_mR*r4 + ηm*P_mR + ηm*Q_mR*ri^2
    M       = ηi*κm - ηi - ηm*κi + ηm
    DEN     = 2.0*ηi^2*κm*r4 - 2.0*ηi^2*κm + ηi*ηm*κi*κm*r4 + ηi*ηm*κi -
                ηi*ηm*κm*r4 + 2.0*ηi*ηm*κm + 2.0*ηi*ηm*r4 - ηi*ηm +
                ηm^2*κi*r4 - ηm^2*κi - ηm^2*r4 + ηm^2
    den1    = ηi*κm*r4 + ηi + ηm*r4 - ηm
    A1_R    = (ηi - ηm)*(r4 - 1.0)*K / DEN
    A1_I    = Q_mI*ri^2*(ηm - ηi)*(r4 - 1.0) / den1
    A1      = A1_R + im*A1_I
    A2_R    = 2.0*(r4 - 1.0)*M*L / (ri^2 * DEN)
    A2_I    = 0.0
    A2      = A2_R + im*A2_I
    A3_R    = (ηi - ηm)*(r4 - 1.0)*(r4 + 1.0)*K / (ri^2 * DEN)
    A3_I    = Q_mI*(ηm - ηi)*(ri^8 - 1.0) / den1
    A3      = A3_R + im*A3_I
    B1_R    = ηi*(κm + 1.0)*L / DEN
    B1_I    = ηi*(κm + 1.0)*(ηi*κm*P_mI*r4 + ηi*P_mI + ηi*Q_mI*ri^2 +
                ηm*P_mI*r4 - ηm*P_mI - ηm*Q_mI*ri^2) / (ηm*(κi + 1.0)*den1)
    B1      = B1_R + im*B1_I
    B2_R    = ηi*ri^2*(κm + 1.0)*K / DEN
    B2_I    = ηi*Q_mI*r4*(κm + 1.0) / den1
    B2      = B2_R + im*B2_I
    λ       = τ + im*σ
    r       = abs(λ)
    #potentials for matrix and inclusion
    if r > ri
        ϕ   =  P_m*(λ + 1.0/λ) + A1/λ
        ϕ′  =  P_m - P_m*λ^-2.0 - A1*λ^-2.0
        ϕ′′ =  2.0*P_m*λ^-3.0 + 2.0*A1*λ^-3.0
        ψ   =  Q_m*(λ + 1.0/λ) + A2/λ + A3*(1.0/(λ^3.0 - λ))
        ψ′  =  Q_m - Q_m*λ^-2.0 - A2*λ^-2.0 - A3*(3.0*λ^2 - 1.0)/((λ^3.0 - λ)^2.0)
        η_loc, κ_loc, ν_loc = ηm, κm, νm
    else
        ϕ   =  B1*(λ + 1.0/λ)
        ϕ′  =  B1 - B1*λ^-2.0
        ϕ′′ =  2.0*B1*λ^-3.0
        ψ   =  B2*(λ + 1.0/λ)
        ψ′  =  B2 - B2*λ^-2.0
        η_loc, κ_loc, ν_loc = ηi, κi, νi
    end
    #field formulas
    ω′      = 1.0 - λ^-2.0
    ω′′     = 2.0*λ^-3.0
    zbar    = conj(λ) + 1.0/conj(λ)          # = conj(z), general point
    Φp      = (ω′*ϕ′′ - ω′′*ϕ′) / ω′^3        # = Φ′(z)
    S       = 4.0*real(ϕ′/ω′)                    # σxx+σyy
    D       = 2.0*( zbar*Φp + ψ′/ω′ )            # σyy-σxx+2iσxy
    sxx     = (S - real(D))/2
    syy     = (S + real(D))/2
    sxy     = imag(D)/2
    szz     = ν_loc*(sxx + syy)
    p       = -2.0*real(1.0/(1.0-λ^-2)*ϕ′)
    conj_ωpr = 1.0 - 1.0/conj(λ)^2         # general conj(ω′(λ))
    vel  = (κ_loc*ϕ - (λ + 1.0/λ)/conj_ωpr*conj(ϕ′) - conj(ψ)) / (2*η_loc)
    return (V   = @SVector([real(vel), imag(vel)]),
            p   = p,
            τ   = @SMatrix([sxx+p  sxy; sxy  syy+p]),
            τzz = szz + p)
end

#Parameter presets
function preset(test)
    if test == :Schmid_ps
        ηm = 1.0; ηi = 1000.0;  ξm, ξi = 1e10, 1e10
        ε̇ = -0.5;  γ̇ = 0.0;  ζ̇ = 0.0;  t = 2.0;  α = 0.0
    elseif test == :Schmid_45deg_ps
        ηm, ηi = 1000, 1.0;  ξm, ξi = 1e10, 1e10
        ε̇ = -0.5;  γ̇ = 0.0;  ζ̇ = 0.0;  t = 2.0;  α = 45/180*pi
    elseif test == :Schmid_30deg_ps
        ηm = 1.0; ηi = 1000.0;  ξm, ξi = 1e10, 1e10
        ε̇ = -0.5;  γ̇ = 0.0;  ζ̇ = 0.0;  t = 2.0;  α = 30/180*pi
    elseif test == :Schmid_ss
        ηm = 1.0; ηi = 1000.0;  ξm, ξi = 1e10, 1e10
        ε̇ = 0.0;  γ̇ = 0.5;  ζ̇ = 0.0;  t = 2.0;  α = 0.0
    elseif test == :Schmid_3
        ηm = 1.0; ηi = 1e-2;  ξm, ξi = 1e10, 1e10
        ε̇ = -1.0;  γ̇ = 0.0;  ζ̇ = 0.0;  t = 1.01;  α = 0.0
    else
        error("unknown test $test")
    end
    return (ηm=ηm, ηi=ηi, ξm=ξm, ξi=ξi, t=t, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇, α=α)
end

#Build solution for plotting
function sample_ring(test; ρ_range, nρ = 400, nθ = 400)
    p    = preset(test)
    ρs   = range(ρ_range[1], ρ_range[2], length = nρ + 1)
    θs   = range(0, 2π,      length = nθ + 1)
    ζ  = [ρ*cis(θ) for θ in θs, ρ in ρs]
    z  = joukowski.(ζ)
    P    = similar(ζ, Float64)
    τmax = similar(ζ, Float64)
    for idx in eachindex(ζ)
        λ = ζ[idx]
        sol = Analytics_new(@SVector([real(λ), imag(λ)]); params = p)
        P[idx]    = sol.p / p.ηm
        sxx, sxy = sol.τ[1,1], sol.τ[1,2]
        syy      = sol.τ[2,2]
        τmax[idx] = sqrt(((sxx - syy)/2)^2 + sxy^2) / p.ηm
    end
    return (; ζ, z, P, τmax)
end

#Plotting
function plot_Joukowsky(test; nρ = 400, nθ = 400, zoom = 4.0)
    p  = preset(test)
    ri = t_to_ri(p.t)
    clast  = sample_ring(test; ρ_range = (1.0, ri),     nρ = nρ, nθ = nθ)
    matrix = sample_ring(test; ρ_range = (ri, zoom*ri), nρ = nρ, nθ = nθ)
    fig = Figure(size = (650, 650))
    ax  = Axis(fig[1, 1]; aspect = DataAspect(), title = "Pressure  p/µₘ  ($(String(test)))")
    hidedecorations!(ax); hidespines!(ax)
    hm = nothing
    for region in (clast, matrix)
        X = real.(region.z)
        Y = imag.(region.z)
        hm = surface!(ax, X, Y, zeros(size(X)); color = region.P, shading = NoShading,
                       colormap = (Reverse(:matter), 1))
    end
    Colorbar(fig[1, 2], hm, label = "p / µₘ")
    save("fig9_style_$(test).png", fig)
    return fig
end

function Analytics_Schmid2003(X; ηm, mc, γ̇, ε̇, α, t)
    ri = t_to_ri(t)          # = rc in the MATLAB script
    τ, σ = X[1], X[2]
    λ = τ + im*σ
    r = abs(λ)
    BC = (2.0*ε̇ - im*γ̇)*exp(2.0*im*α)
    ReB, ImB = real(BC), imag(BC)
    B1 = ri^4*mc + ri^4 - 1.0 + mc
    B2 = ri^4*mc + ri^4 - mc + 1.0
    B3 = ri^4*mc - mc - ri^4 + 1.0
    B4 = -ri^4*mc - mc - ri^4 + 1.0
    B5 = ri^8*mc - mc - ri^8 + 1.0
    C1 = B3*ri^2*(im*ImB/B1 - ReB/B2)     # coeff of 1/ζ    in φ_m
    C2 = B5*(im*ImB/B1 - ReB/B2)          # coeff of 1/(ζ³-ζ) in ψ_m
    Bc1 = im*mc*B4*γ̇/(2.0*B1) - ri^2*(mc-1.0)*(im*mc*ImB/B1 - ReB/B2)   # coeff of ω(ζ) in φ_c
    Bc2 = -2.0*mc*ri^4*(im*ImB/B1 + ReB/B2)                              # coeff of ω(ζ) in ψ_c
    if r > ri
        ϕ   = -im/2*γ̇*(λ + 1.0/λ) + C1/λ
        ϕ′  = -im/2*γ̇*(1.0 - λ^-2.0) - C1*λ^-2.0
        ϕ′′ = -im/2*γ̇*(2.0*λ^-3.0) + 2.0*C1*λ^-3.0
        ψ   = -(ReB + im*ImB)*(λ + 1.0/λ) + C2*(1.0/(λ^3.0 - λ))
        ψ′  = -(ReB + im*ImB)*(1.0 - λ^-2.0) - C2*(3.0*λ^2 - 1.0)/((λ^3.0 - λ)^2.0)
        η_loc = ηm
    else
        ϕ   = Bc1*(λ + 1.0/λ)
        ϕ′  = Bc1*(1.0 - λ^-2.0)
        ϕ′′ = Bc1*(2.0*λ^-3.0)
        ψ   = Bc2*(λ + 1.0/λ)
        ψ′  = Bc2*(1.0 - λ^-2.0)
        η_loc = ηm*mc
    end
    #field solutions
    ω′  = 1.0 - λ^-2.0
    ω′′ = 2.0*λ^-3.0
    zbar = conj(λ) + 1.0/conj(λ)
    Φp   = (ω′*ϕ′′ - ω′′*ϕ′) / ω′^3
    S = 4.0*real(ϕ′/ω′)
    D = 2.0*( zbar*Φp + ψ′/ω′ )
    sxx = (S - real(D))/2
    syy = (S + real(D))/2
    sxy = imag(D)/2
    p   = -2.0*real(1.0/(1.0 - λ^-2.0)*ϕ′) * ηm   
    conj_ωpr = 1.0 - 1.0/conj(λ)^2
    vel = (ϕ - (λ + 1.0/λ)/conj_ωpr*conj(ϕ′) - conj(ψ)) / (2.0*η_loc) * ηm
    return (V   = @SVector([real(vel), imag(vel)]),
            p   = p,
            τ   = @SMatrix([sxx*ηm+p  sxy*ηm; sxy*ηm  syy*ηm+p]),
            τmax = sqrt(((sxx-syy)/2)^2 + sxy^2) * ηm)
end

function sample_ring_compare(test; ρ_range, nρ = 300, nθ = 300)
    p  = preset(test)
    ri = t_to_ri(p.t)
    α_schmid = -p.α
    ρs = range(ρ_range[1], ρ_range[2], length = nρ + 1)
    θs = range(0, 2π,      length = nθ + 1)
    ζ  = [ρ*cis(θ) for θ in θs, ρ in ρs]
    z  = joukowski.(ζ)
    P_mine, Vx_mine, Vy_mine = similar(ζ, Float64), similar(ζ, Float64), similar(ζ, Float64)
    P_ref,  Vx_ref,  Vy_ref  = similar(ζ, Float64), similar(ζ, Float64), similar(ζ, Float64)

    for idx in eachindex(ζ)
        λ = ζ[idx]
        X = @SVector([real(λ), imag(λ)])
        sol_mine = Analytics_new(X; params = p)
        P_mine[idx]  = sol_mine.p / p.ηm
        Vx_mine[idx] = sol_mine.V[1]
        Vy_mine[idx] = sol_mine.V[2]
        sol_ref = Analytics_Schmid2003(X; ηm = p.ηm, mc = p.ηi/p.ηm,
                                        γ̇ = p.γ̇, ε̇ = p.ε̇, α = α_schmid, t = p.t)
        P_ref[idx]  = sol_ref.p / p.ηm
        Vx_ref[idx] = sol_ref.V[1]
        Vy_ref[idx] = sol_ref.V[2]
    end

    return (; z, P_mine, Vx_mine, Vy_mine, P_ref, Vx_ref, Vy_ref)
end

function plot_comparison(test; nρ = 300, nθ = 300, zoom = 4.0)
    p  = preset(test)
    ri = t_to_ri(p.t)

    clast  = sample_ring_compare(test; ρ_range = (1.0, ri),     nρ = nρ, nθ = nθ)
    matrix = sample_ring_compare(test; ρ_range = (ri, zoom*ri), nρ = nρ, nθ = nθ)

    fig = Figure(size = (1500, 1400))
    Label(fig[0, 1:3], "Comparison: Analytics_new vs Analytics_Schmid2003  ($(String(test)))";
          fontsize = 20, font = :bold)

    rows = [
        ("Vx",  r -> r.Vx_mine, r -> r.Vx_ref),
        ("Vy",  r -> r.Vy_mine, r -> r.Vy_ref),
        ("P",   r -> r.P_mine,  r -> r.P_ref),
    ]

    for (i, (label, mine_fn, ref_fn)) in enumerate(rows)
        clast_mine, clast_ref   = mine_fn(clast),  ref_fn(clast)
        matrix_mine, matrix_ref = mine_fn(matrix), ref_fn(matrix)

        clast_diff  = clast_mine  .- clast_ref
        matrix_diff = matrix_mine .- matrix_ref

        panels = [
            ("$(label) -- DM_Joukowsky",       clast_mine,  matrix_mine),
            ("$(label) -- Schmid_Joukowsky",     clast_ref,   matrix_ref),
            ("$(label) -- Difference", clast_diff,  matrix_diff),
        ]

        for (j, (ttl, cvals, mvals)) in enumerate(panels)
            gl = GridLayout(fig[i, j])
            ax = Axis(gl[1, 1]; aspect = DataAspect(), title = ttl)
            hidedecorations!(ax); hidespines!(ax)

            hm = nothing
            for (region, vals) in ((clast, cvals), (matrix, mvals))
                X = real.(region.z); Y = imag.(region.z)
                hm = surface!(ax, X, Y, zeros(size(X)); color = vals, shading = NoShading,
                               colormap = (Reverse(:matter), 1))
            end
            Colorbar(gl[2, 1], hm; vertical = false, flipaxis = false, height = 12)
            rowsize!(gl, 2, Relative(0.06))
        end
    end

    save("comparison_Joukowsky_$(test).png", fig)
    return fig
end

plot_comparison(:Schmid_ps)

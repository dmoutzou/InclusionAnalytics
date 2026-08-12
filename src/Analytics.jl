# # Analytical field at a NORMALISED point Xn ∈ the unit box.
# function f_anal(Xn, params; geometry=:circular)
#     if geometry == :circular
#         return analytics_circle
#     elseif geometry == :elliptical
#         return analytics_ellipse(Xn; params=params)
#     else
#         error("Unknown geometry: $geometry")
#     end
# end

# New analytical solution (including OOP)
function analytics_circle(X; params)
    ηm, ηi, ξm, ξi, ri, r2, t, α, sc, ε̇, γ̇, ζ̇, ε̇zz = params
    x, y = X[1], X[2]
    νm = (3ξm - 2ηm) / (2*(3ξm + ηm))
    νi = (3ξi - 2ηi) / (2*(3ξi + ηi))
    κm = 3 - 4νm
    κi = 3 - 4νi
    Em = 2ηm*(1 + νm)
    Ei = 2ηi*(1 + νi)

    P_mR = 2*ηm*(ζ̇ + νm*ε̇zz) / (κm - 1)
    P_mI = -ηm*γ̇ / (κm + 1)
    P_m  = P_mR + im*P_mI
    Q_m  = (im*γ̇ - 2ε̇)*ηm

    D   = 2ηi + (κi - 1)*ηm
    B1R = ((κm + 1)*ηi*P_mR + 2*ηi*ηm*(νi - νm)*ε̇zz) / D
    B1I = (κm + 1)*ηi*P_mI / ((κi + 1)*ηm)
    B1  = B1R + im*B1I

    A1  = (ηi - ηm)*conj(Q_m)*ri^2 / (κm*ηi + ηm)
    A3  = ri^2 * A1
    B2  = conj(A1/ri^2 + conj(Q_m))
    A2  = 2*ri^2 * (B1R - P_mR)

    z   = x + im*y
    z̄   = conj(z)
    r   = abs(z)
    if r > ri
        ϕ   =  P_m*z  + A1/z
        ϕ′  =  P_m    - A1/z^2
        ϕ′′ =           2A1/z^3
        ψ   =  Q_m*z  + A2/z   + A3/z^3
        ψ′  =  Q_m    - A2/z^2 - 3A3/z^4
        η_loc, κ_loc, ν_loc, E_loc = ηm, κm, νm, Em
    else
        ϕ   =  B1*z
        ϕ′  =  B1  + 0im
        ϕ′′ =  0.0 + 0im
        ψ   =  B2*z
        ψ′  =  B2  + 0im
        η_loc, κ_loc, ν_loc, E_loc = ηi, κi, νi, Ei
    end
    Qval = z̄*ϕ′′ + ψ′
    sxx  = 2real(ϕ′) - real(Qval)
    syy  = 2real(ϕ′) + real(Qval)
    sxy  = imag(Qval)
    szz  = ν_loc*(sxx + syy) + E_loc*ε̇zz
    p    = -(sxx + syy + szz) / 3
    vel  = (κ_loc*ϕ - z*conj(ϕ′) - conj(ψ)) / (2*η_loc) - ν_loc*ε̇zz*z
    return (V   = @SVector([real(vel), imag(vel)]),
            p   = p,
            τ   = @SMatrix([sxx+p  sxy; sxy  syy+p]),
            τzz = szz + p,
            σ   = @SMatrix([sxx  sxy; sxy  syy]),
            σzz = szz)
end

# physical position from the conformal coordinate
joukowski(λ) = λ + 1.0/λ          # ω(λ) = z   (change here if your map differs)
t_to_ri(t) = sqrt((t - 1.0)*(t + 1.0)) / (t - 1.0)

# The numerical solver is using the true cartesian coordinates, so I need the inverse mapping of the Joukowsky transform
function inv_joukowski(z)
    disc = sqrt(z^2 - 4.0 + 0im)
    ζ1   = (z + disc)/2
    ζ2   = (z - disc)/2
    return abs(ζ1) >= abs(ζ2) ? ζ1 : ζ2
end
to_zeta(X) = (ζ = inv_joukowski(complex(X[1], X[2])); @SVector([real(ζ), imag(ζ)]))

# Function to acquire the axes of the inclusion based on t (
function ellipse_axes(t)
    ri = t_to_ri(t)
    a  = ri + 1.0/ri
    b  = ri - 1.0/ri
    return a, b
end

#New analytical solution for ellipse (including OOP)
function analytics_ellipse(X; params)
    ηm, ηi, ξm, ξi, ri, r2, t, α, sc, ε̇, γ̇, ζ̇, ε̇zz = params
    ri      = t_to_ri(t)
    Ζ       = to_zeta(sc.*X)
    τ, σ    = Ζ[1], Ζ[2]
    # Kolosov
    νm      = (3.0ξm - 2.0ηm) / (2.0*(3.0ξm + ηm))
    νi      = (3.0ξi - 2.0ηi) / (2.0*(3.0ξi + ηi))
    κm      = 3 - 4νm
    κi      = 3 - 4νi
    Em      = 2.0*ηm*(1.0 + νm)
    Ei      = 2.0*ηi*(1.0 + νi)
    # Boundary terms with angle of rotation
    P_mR    = 2.0*ηm*(ζ̇ + νm*ε̇zz) / (κm - 1.0)
    P_mI    = -ηm*γ̇ / (κm + 1.0)
    P_m     = P_mR + im*P_mI
    Q_m     = (im*γ̇ - 2.0ε̇)*ηm*exp(-2.0*im*α)
    Q_mR    = real(Q_m); Q_mI = imag(Q_m)
    # radius factor
    r4      = ri^4
    # out-of-plane forcing of the velocity-continuity condition
    Δzz     = 2.0*ηi*ηm*ε̇zz*(νm - νi)
    # Solution coefficients
    K       = -2.0*ηi*κm*P_mR + 2.0*ηi*P_mR + 2.0*ηi*Q_mR*ri^2 + 2.0*ηm*κi*P_mR +
                ηm*κi*Q_mR*ri^2 - 2.0*ηm*P_mR - ηm*Q_mR*ri^2 + 2.0*Δzz
    L       = ηi*κm*P_mR*r4 - ηi*P_mR - ηi*Q_mR*ri^2 + ηm*P_mR*r4 + ηm*P_mR + ηm*Q_mR*ri^2
    M       = ηi*κm - ηi - ηm*κi + ηm
    DEN     = 2.0*ηi^2*κm*r4 - 2.0*ηi^2*κm + ηi*ηm*κi*κm*r4 + ηi*ηm*κi -
                ηi*ηm*κm*r4 + 2.0*ηi*ηm*κm + 2.0*ηi*ηm*r4 - ηi*ηm +
                ηm^2*κi*r4 - ηm^2*κi - ηm^2*r4 + ηm^2
    den1    = ηi*κm*r4 + ηi + ηm*r4 - ηm
    den2    = ηi*κm*r4 - ηi + ηm*r4 + ηm
    A1_R    = (ηi - ηm)*(r4 - 1.0)*K / DEN
    A1_I    = Q_mI*ri^2*(ηm - ηi)*(r4 - 1.0) / den1
    A1      = A1_R + im*A1_I
    A2_R    = 2.0*(r4 - 1.0)*(M*L - Δzz*den2) / (ri^2 * DEN)
    A2_I    = 0.0
    A2      = A2_R + im*A2_I
    A3_R    = (ηi - ηm)*(r4 - 1.0)*(r4 + 1.0)*K / (ri^2 * DEN)
    A3_I    = Q_mI*(ηm - ηi)*(ri^8 - 1.0) / den1
    A3      = A3_R + im*A3_I
    B1_R    = (ηi*(κm + 1.0)*L - Δzz*den1) / DEN
    B1_I    = ηi*(κm + 1.0)*(ηi*κm*P_mI*r4 + ηi*P_mI + ηi*Q_mI*ri^2 +
                ηm*P_mI*r4 - ηm*P_mI - ηm*Q_mI*ri^2) / (ηm*(κi + 1.0)*den1)
    B1      = B1_R + im*B1_I
    B2_R    = ηi*ri^2*(κm + 1.0)*K / DEN
    B2_I    = ηi*Q_mI*r4*(κm + 1.0) / den1
    B2      = B2_R + im*B2_I
    λ       = τ + im*σ
    r       = abs(λ)
    # potentials for matrix and inclusion
    if r > ri
        ϕ   =  P_m*(λ + 1.0/λ) + A1/λ
        ϕ′  =  P_m - P_m*λ^-2.0 - A1*λ^-2.0
        ϕ′′ =  2.0*P_m*λ^-3.0 + 2.0*A1*λ^-3.0
        ψ   =  Q_m*(λ + 1.0/λ) + A2/λ + A3*(1.0/(λ^3.0 - λ))
        ψ′  =  Q_m - Q_m*λ^-2.0 - A2*λ^-2.0 - A3*(3.0*λ^2 - 1.0)/((λ^3.0 - λ)^2.0)
        η_loc, κ_loc, ν_loc, E_loc = ηm, κm, νm, Em
    else
        ϕ   =  B1*(λ + 1.0/λ)
        ϕ′  =  B1 - B1*λ^-2.0
        ϕ′′ =  2.0*B1*λ^-3.0
        ψ   =  B2*(λ + 1.0/λ)
        ψ′  =  B2 - B2*λ^-2.0
        η_loc, κ_loc, ν_loc, E_loc = ηi, κi, νi, Ei
    end
    # field formulas
    ω        = λ + 1.0/λ                       # = z (scaled)
    ω′       = 1.0 - λ^-2.0
    ω′′      = 2.0*λ^-3.0
    zbar     = conj(λ) + 1.0/conj(λ)           # = conj(z)
    Φp       = (ω′*ϕ′′ - ω′′*ϕ′) / ω′^3        # = Φ′(z)
    S        = 4.0*real(ϕ′/ω′)                 # σxx+σyy
    D        = 2.0*( zbar*Φp + ψ′/ω′ )         # σyy-σxx+2iσxy
    sxx      = (S - real(D))/2
    syy      = (S + real(D))/2
    sxy      = imag(D)/2
    szz      = ν_loc*(sxx + syy) + E_loc*ε̇zz
    p        = -(sxx + syy + szz) / 3.0
    conj_ωpr = 1.0 - 1.0/conj(λ)^2
    vel      = (κ_loc*ϕ - ω/conj_ωpr*conj(ϕ′) - conj(ψ)) / (2*η_loc) - ν_loc*ε̇zz*ω
    return (V   = @SVector([real(vel)/sc, imag(vel)/sc]),
            p   = p,
            τ   = @SMatrix([sxx+p  sxy; sxy  syy+p]),
            τzz = szz + p,
            σ   = @SMatrix([sxx  sxy; sxy  syy]),
            σzz = szz)
end

#Analytical solution for ellipse from Schmid2003
function Analytics_Schmid2003(X; ηm, mc, γ̇, ε̇, α, t)
    ri = t_to_ri(t)
    τ, σ = X[1], X[2]
    λ = τ + im*σ
    r = abs(λ)
    BC = (2.0*ε̇ - im*γ̇)*exp(-2.0*im*α)
    ReB, ImB = real(BC), imag(BC)
    B1 = ri^4*mc + ri^4 - 1.0 + mc
    B2 = ri^4*mc + ri^4 - mc + 1.0
    B3 = ri^4*mc - mc - ri^4 + 1.0
    B4 = -ri^4*mc - mc - ri^4 + 1.0
    B5 = ri^8*mc - mc - ri^8 + 1.0
    C1 = B3*ri^2*(im*ImB/B1 - ReB/B2)
    C2 = B5*(im*ImB/B1 - ReB/B2)
    Bc1 = im*mc*B4*γ̇/(2.0*B1) - ri^2*(mc-1.0)*(im*mc*ImB/B1 - ReB/B2)
    Bc2 = -2.0*mc*ri^4*(im*ImB/B1 + ReB/B2)
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
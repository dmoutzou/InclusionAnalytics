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

function analytics_ellipse(X; params)
    ηm, ηi, ξm, ξi, t, γ̇, ε̇, ζ̇ , α = params
    ri      = t_to_ri(t)
    τ, σ    = X[1], X[2]
    # Kolosov
    νm      = (3.0ξm - 2.0ηm) / (2.0*(3.0ξm + ηm))
    νi      = (3.0ξi - 2.0ηi) / (2.0*(3.0ξi + ηi))
    κm      = 3 - 4νm
    κi      = 3 - 4νi
    # Boundary terms with angle of rotation
    P_m     = -im*ηm*γ̇/(κm+1.0) + 2.0*ηm*ζ̇/(κm-1)
    Q_m     = (im*γ̇ - 2.0ε̇)*ηm*exp(-2.0*im*α)
    P_mR    = real(P_m); P_mI = imag(P_m)
    Q_mR    = real(Q_m); Q_mI = imag(Q_m)
    # radius factor
    r4      = ri^4
    # Solution coefficients
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
    # potentials for matrix and inclusion
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
    # field formulas
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
    p       = -(sxx + syy + szz) / 3.0        # full 3D mean pressure -(σxx+σyy+σzz)/3
    conj_ωpr = 1.0 - 1.0/conj(λ)^2         # general conj(ω′(λ))
    vel  = (κ_loc*ϕ - (λ + 1.0/λ)/conj_ωpr*conj(ϕ′) - conj(ψ)) / (2*η_loc)
    return (V   = @SVector([real(vel), imag(vel)]),
            p   = p,
            τ   = @SMatrix([sxx+p  sxy; sxy  syy+p]),
            τzz = szz + p)
end
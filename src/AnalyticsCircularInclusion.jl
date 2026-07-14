# Old analytical solution (Duretz)
# Matrix potentials
ϕ_mat(z,   k, γ̇, ε̇, ζ̇, ηm, κm, ri) =  2*ηm*ζ̇/(κm-1)*z - im/2*ηm*γ̇*z      -   (im*γ̇ + 2*ε̇) * ηm * k[1] * ri^2 / z^1
ϕ′_mat(z,  k, γ̇, ε̇, ζ̇, ηm, κm, ri) =  2*ηm*ζ̇/(κm-1)   - im/2*ηm*γ̇        +   (im*γ̇ + 2*ε̇) * ηm * k[1] * ri^2 / z^2
ϕ′′_mat(z, k, γ̇, ε̇, ζ̇, ηm, κm, ri) =                  - 2*(im*γ̇ + 2*ε̇) * ηm * k[1] * ri^2 / z^3
ψ_mat(z,   k, γ̇, ε̇, ζ̇, ηm, κm, ri) =                   (im*γ̇ - 2*ε̇)*ηm*z -   (im*γ̇ + 2*ε̇) * ηm * k[1] * ri^4 / z^3
ψ′_mat(z,  k, γ̇, ε̇, ζ̇, ηm, κm, ri) =                   (im*γ̇ - 2*ε̇)*ηm   + 3*(im*γ̇ + 2*ε̇) * ηm * k[1] * ri^4 / z^4
# Inclusion potentials
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
    k1 = (ηi - ηm) / (ηm + ηi*κm)
    k2 = ηm*(κm + 1) / (2*(ηm + ηi*κm))
    k3 = (Mm - Mi) / (Mi + ηm)
    k4 = 1 + k3
    return @SVector([k1, k2, k3, k4])
end

function Analytics_old(X; params)
    ηm, ηi, ξm, ξi, γ̇, ε̇, ζ̇, ri = params
    z   = X[1] + im*X[2]
    z̄   = conj(z)
    r   = abs(z)
    θ   = angle(z)
    νm  = (3*ξm - 2*ηm)/(2*(3*ξm + ηm))
    νi  = (3*ξi - 2*ηi)/(2*(3*ξi + ηi))
    κm  = 3 - 4*νm
    κi  = 3 - 4*νi
    k   = constants_analytical(ri, ηm, ηi, ξm, ξi, γ̇, ε̇, ζ̇)

    if r > ri
        η, κ, ν = ηm, κm, νm
        ϕ′  = ϕ′_mat(z,  k, γ̇, ε̇, ζ̇, ηm, κm, ri)
        ϕ′′ = ϕ′′_mat(z, k, γ̇, ε̇, ζ̇, ηm, κm, ri)
        ϕ   = ϕ_mat(z,   k, γ̇, ε̇, ζ̇, ηm, κm, ri)
        ψ   = ψ_mat(z,   k, γ̇, ε̇, ζ̇, ηm, κm, ri)
        ψ′  = ψ′_mat(z,  k, γ̇, ε̇, ζ̇, ηm, κm, ri)
        W_pert   =  k[3] * ζ̇ * ri^2/z̄
        fac      =  2 * ηm * k[3] * ζ̇ * ri^2/r^2
        sxx_pert = -fac*cos(2θ)
        syy_pert =  fac*cos(2θ)
        sxy_pert = -fac*sin(2θ)
    else
        η, κ, ν = ηi, κi, νi
        ϕ′  = ϕ′_inc(z,  k, γ̇, ε̇, ζ̇, ηi, κi, ri)
        ϕ′′ = ϕ′′_inc(z, k, γ̇, ε̇, ζ̇, ηi, κi, ri)
        ϕ   = ϕ_inc(z,   k, γ̇, ε̇, ζ̇, ηi, κi, ri)
        ψ   = ψ_inc(z,   k, γ̇, ε̇, ζ̇, ηi, κi, ri)
        ψ′  = ψ′_inc(z,  k, γ̇, ε̇, ζ̇, ηi, κi, ri)
        W_pert   = 0.0 + 0.0im
        sxx_pert = 0.0; syy_pert = 0.0; sxy_pert = 0.0
    end

    Q    = z̄*ϕ′′ + ψ′
    sxx  = 2*real(ϕ′) - real(Q) + sxx_pert
    syy  = 2*real(ϕ′) + real(Q) + syy_pert
    sxy  = imag(Q)              + sxy_pert
    szz  = ν*(sxx + syy)
    p    = -1/3*(sxx + syy + szz)
    W    = 1/(2*η)*(κ*ϕ - z*conj(ϕ′) - conj(ψ)) + W_pert
    return (V   = @SVector([real(W), imag(W)]),
            p   = p,
            τ   = @SMatrix([sxx+p  sxy; sxy  syy+p]),
            τzz = szz + p)
end

# New analytical solution (Moutzouris_2026)
function Analytics_new(X; params)
    ηm, ηi, ξm, ξi, γ̇, ε̇, ζ̇, ri = params
    x, y = X[1], X[2]
    νm = (3ξm - 2ηm) / (2*(3ξm + ηm))
    νi = (3ξi - 2ηi) / (2*(3ξi + ηi))
    κm = 3 - 4νm
    κi = 3 - 4νi
    P_m  = -im*ηm*γ̇/(κm+1) + 2*ηm*ζ̇/(κm-1)
    Q_m  = (im*γ̇ - 2ε̇)*ηm
    P_mR = real(P_m)
    P_mI = imag(P_m)
    B1R = (κm+1)*ηi*P_mR / (2ηi + (κi-1)*ηm)
    B1I = (κm+1)*ηi*P_mI / ((κi+1)*ηm)
    B1  = B1R + im*B1I
    A1  = (ηi - ηm)*conj(Q_m)*ri^2 / (κm*ηi + ηm)
    A3  = ri^2 * A1
    B2  = conj(A1/ri^2 + conj(Q_m))
    A2  = 2ri^2 * ((κm-1)*ηi - (κi-1)*ηm) / (2ηi + (κi-1)*ηm) * P_mR
    z   = x + im*y
    z̄   = conj(z)
    r   = abs(z)
    if r > ri
        ϕ   =  P_m*z  + A1/z
        ϕ′  =  P_m    - A1/z^2
        ϕ′′ =           2A1/z^3
        ψ   =  Q_m*z  + A2/z   + A3/z^3
        ψ′  =  Q_m    - A2/z^2 - 3A3/z^4
        η_loc, κ_loc, ν_loc = ηm, κm, νm
    else
        ϕ   =  B1*z
        ϕ′  =  B1  + 0im
        ϕ′′ =  0.0 + 0im
        ψ   =  B2*z
        ψ′  =  B2  + 0im
        η_loc, κ_loc, ν_loc = ηi, κi, νi
    end
    Qval = z̄*ϕ′′ + ψ′
    sxx  = 2real(ϕ′) - real(Qval)
    syy  = 2real(ϕ′) + real(Qval)
    sxy  = imag(Qval)
    szz  = ν_loc*(sxx + syy)
    p    = -(sxx + syy + szz) / 3
    vel  = (κ_loc*ϕ - z*conj(ϕ′) - conj(ψ)) / (2*η_loc)
    return (V   = @SVector([real(vel), imag(vel)]),
            p   = p,
            τ   = @SMatrix([sxx+p  sxy; sxy  syy+p]),
            τzz = szz + p)
end
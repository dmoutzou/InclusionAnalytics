function Analytics_new(X; params)
    ηm, ηi, ξm, ξi, ri, γ̇, ε̇, ζ̇ = params
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

function Analytical_Solution_new(x, y, test)
    if test == :Schmid2003
        ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e10, 1e10
        ε̇ =-1e0;  γ̇ = 0.0;  ζ̇ = 0.0;  ri = 0.1
    elseif test == :Duretz2026_1
        ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e0, 1e0
        ε̇ =-1e0;  γ̇ = 0.0;  ζ̇ = 0.0;  ri = 0.1
    elseif test == :Duretz2026_2
        ηm, ηi = 1.0, 1e-1;  ξm, ξi = 1e0, 1e0
        ε̇ = 1e-10;  γ̇ = 0.0;  ζ̇ = 1e0;  ri = 0.1
    elseif test == :Duretz2026_3
        ηm, ηi = 1.0, 1e-1;  ξm, ξi = 1e2, 1e2
        ε̇ = 1e-10;  γ̇ = 1.0;  ζ̇ = 1e0;  ri = 0.1
    end
    params = (ηm=ηm, ηi=ηi, ξm=ξm, ξi=ξi, ri=ri, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇)
    sol = Analytics_new(@SVector([x, y]); params=params)
    return sol.p, sol.V[1], sol.V[2]
end

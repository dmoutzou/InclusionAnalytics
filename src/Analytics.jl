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

# Analytical solution for ellipse from Schmid2003  --> EFS
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
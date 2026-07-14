using CairoMakie, StaticArrays, LinearAlgebra

let

# ============================================================
# Parameters (Duretz2026_2 test case, shared by both files)
# ============================================================
ηm = 1.0;  ηi = 0.1
ξm = 1.0;  ξi = 1.0
ri = 0.1;  γ̇ = 1.0;  ε̇ = 1e-10;  ζ̇ = 1.0

νm = (3ξm - 2ηm) / (2*(3ξm + ηm))
νi = (3ξi - 2ηi) / (2*(3ξi + ηi))
κm = 3 - 4νm
κi = 3 - 4νi
Mm = ξm + ηm/3
Mi = ξi + ηi/3

# ============================================================
# OLD formulation — Analytics_v12.jl
# ============================================================
# Deviatoric & volumetric constants
k1 = (ηi - ηm) / (ηm + ηi*κm)
k2 = ηm*(κm + 1) / (2*(ηm + ηi*κm))
k3 = (Mm - Mi) / (Mi + ηm)
k4 = 1 + k3
k  = (k1, k2, k3, k4)

ϕ_mat(z)   =  2ηm*ζ̇/(κm-1)*z - im/2*ηm*γ̇*z      - (im*γ̇ + 2ε̇)*ηm*k1*ri^2/z
ϕ′_mat(z)  =  2ηm*ζ̇/(κm-1)   - im/2*ηm*γ̇         + (im*γ̇ + 2ε̇)*ηm*k1*ri^2/z^2
ϕ′′_mat(z) =                   - 2*(im*γ̇ + 2ε̇)*ηm*k1*ri^2/z^3
ψ_mat(z)   =                   (im*γ̇ - 2ε̇)*ηm*z   - (im*γ̇ + 2ε̇)*ηm*k1*ri^4/z^3
ψ′_mat(z)  =                   (im*γ̇ - 2ε̇)*ηm     + 3*(im*γ̇ + 2ε̇)*ηm*k1*ri^4/z^4

ϕ_inc(z)   =  2ηi*k4*ζ̇/(κi-1)*z - im/2*ηi*γ̇*z
ϕ′_inc(z)  =  2ηi*k4*ζ̇/(κi-1)   - im/2*ηi*γ̇
ϕ′′_inc(z) =  0.0 + 0im
ψ_inc(z)   =  2*(im*γ̇ - 2ε̇)*ηi*k2*z
ψ′_inc(z)  =  2*(im*γ̇ - 2ε̇)*ηi*k2

function old_solution(x, y)
    z = x + im*y
    z̄ = conj(z)
    r = abs(z)
    θ = angle(z)
    if r > ri
        η, κ, ν = ηm, κm, νm
        ϕ, ϕ′, ϕ′′, ψ, ψ′ = ϕ_mat(z), ϕ′_mat(z), ϕ′′_mat(z), ψ_mat(z), ψ′_mat(z)
        fac      =  2ηm*k3*ζ̇*ri^2/r^2
        sxx_pert = -fac*cos(2θ)
        syy_pert =  fac*cos(2θ)
        sxy_pert = -fac*sin(2θ)
        W_pert   =  k3*ζ̇*ri^2/z̄
    else
        η, κ, ν = ηi, κi, νi
        ϕ, ϕ′, ϕ′′, ψ, ψ′ = ϕ_inc(z), ϕ′_inc(z), ϕ′′_inc(z), ψ_inc(z), ψ′_inc(z)
        sxx_pert = syy_pert = sxy_pert = 0.0
        W_pert   = 0.0 + 0im
    end
    Q    = z̄*ϕ′′ + ψ′
    sxx  = 2real(ϕ′) - real(Q) + sxx_pert
    syy  = 2real(ϕ′) + real(Q) + syy_pert
    sxy  = imag(Q)              + sxy_pert
    szz  = ν*(sxx + syy)
    p    = -(sxx + syy + szz)/3
    W    = (κ*ϕ - z*conj(ϕ′) - conj(ψ))/(2η) + W_pert
    return p, real(W), imag(W), sxx, syy, sxy
end

# ============================================================
# NEW formulation — New_formulation_v3.jl
# ============================================================
P_m  = -im*ηm*γ̇/(κm+1) + 2ηm*ζ̇/(κm-1)
Q_m  = (im*γ̇ - 2ε̇)*ηm
P_mR = real(P_m);  P_mI = imag(P_m)

B1R  = (κm+1)*ηi*P_mR / (2ηi + (κi-1)*ηm)
B1I  = (κm+1)*ηi*P_mI / ((κi+1)*ηm)
B1   = B1R + im*B1I
A1   = (ηi - ηm)*conj(Q_m)*ri^2 / (κm*ηi + ηm)
A3   = ri^2 * A1
B2   = conj(A1/ri^2 + conj(Q_m))
A2   = 2ri^2 * ((κm-1)*ηi - (κi-1)*ηm) / (2ηi + (κi-1)*ηm) * P_mR

function new_solution(x, y)
    z = x + im*y
    z̄ = conj(z)
    r = abs(z)
    if r > ri
        ϕ   = P_m*z + A1/z
        ϕ′  = P_m   - A1/z^2
        ϕ′′ =         2A1/z^3
        ψ   = Q_m*z + A2/z  + A3/z^3
        ψ′  = Q_m   - A2/z^2 - 3A3/z^4
        η_loc, κ_loc, ν_loc = ηm, κm, νm
    else
        ϕ   = B1*z
        ϕ′  = B1   + 0im
        ϕ′′ = 0.0  + 0im
        ψ   = B2*z
        ψ′  = B2   + 0im
        η_loc, κ_loc, ν_loc = ηi, κi, νi
    end
    Qval = z̄*ϕ′′ + ψ′
    sxx  = 2real(ϕ′) - real(Qval)
    syy  = 2real(ϕ′) + real(Qval)
    sxy  = imag(Qval)
    szz  = ν_loc*(sxx + syy)
    p    = -(sxx + syy + szz)/3
    vel  = (κ_loc*ϕ - z*conj(ϕ′) - conj(ψ))/(2η_loc)
    return p, real(vel), imag(vel), sxx, syy, sxy
end

# ============================================================
# Grid evaluation
# ============================================================
n  = 200
x  = LinRange(-0.3, 0.3, n)
y  = LinRange(-0.3, 0.3, n)

fields = [:P, :Vx, :Vy, :Sxx, :Syy, :Sxy]
OLD = NamedTuple{Tuple(fields)}(zeros(n,n) for _ in fields)
NEW = NamedTuple{Tuple(fields)}(zeros(n,n) for _ in fields)

for (j, yj) in enumerate(y), (i, xi) in enumerate(x)
    p, vx, vy, sxx, syy, sxy = old_solution(xi, yj)
    OLD.P[i,j]   = p;  OLD.Vx[i,j] = vx;  OLD.Vy[i,j] = vy
    OLD.Sxx[i,j] = sxx; OLD.Syy[i,j] = syy; OLD.Sxy[i,j] = sxy

    p, vx, vy, sxx, syy, sxy = new_solution(xi, yj)
    NEW.P[i,j]   = p;  NEW.Vx[i,j] = vx;  NEW.Vy[i,j] = vy
    NEW.Sxx[i,j] = sxx; NEW.Syy[i,j] = syy; NEW.Sxy[i,j] = sxy
end

# ============================================================
# Diagnostics
# ============================================================
println("=== Max |difference| between OLD and NEW ===")
for f in fields
    d = maximum(abs.(OLD[f] .- NEW[f]))
    println("  $f : max|Δ| = $d")
end

# ============================================================
# Plot: OLD | NEW | DIFFERENCE  for P, Vx, Vy
# ============================================================
labels   = ["Pressure  p", "Velocity  Vx", "Velocity  Vy", "Stress  σxx", "Stress  σyy", "Stress  σxy"]
old_data = [OLD.P,  OLD.Vx,  OLD.Vy,  OLD.Sxx,  OLD.Syy,  OLD.Sxy ]
new_data = [NEW.P,  NEW.Vx,  NEW.Vy,  NEW.Sxx,  NEW.Syy,  NEW.Sxy ]

fig = Figure(size=(1300, 1050), fontsize=12)
Label(fig[0, 1:2], "Analytics_v12  (old)", tellwidth=false, fontsize=14, font=:bold)
Label(fig[0, 3:4], "New_formulation_v3  (new)", tellwidth=false, fontsize=14, font=:bold)
Label(fig[0, 5:6], "Difference  (old − new)", tellwidth=false, fontsize=14, font=:bold)

for (row, (lbl, od, nd)) in enumerate(zip(labels, old_data, new_data))
    diff = od .- nd
    clim_sym = maximum(abs, od)
    clim_sym = clim_sym == 0 ? 1.0 : clim_sym
    clim_diff = max(maximum(abs, diff), 1e-12 * clim_sym)

    ax1 = Axis(fig[row, 1], aspect=DataAspect(), title=(row==1 ? "Old" : ""), ylabel=lbl)
    hm1 = heatmap!(ax1, x, y, od; colormap=Reverse(:RdBu), colorrange=(-clim_sym, clim_sym))
    Colorbar(fig[row, 2], hm1, width=12)

    ax2 = Axis(fig[row, 3], aspect=DataAspect(), title=(row==1 ? "New" : ""))
    hm2 = heatmap!(ax2, x, y, nd; colormap=Reverse(:RdBu), colorrange=(-clim_sym, clim_sym))
    Colorbar(fig[row, 4], hm2, width=12)

    ax3 = Axis(fig[row, 5], aspect=DataAspect(), title=(row==1 ? "Difference (old − new)" : ""))
    hm3 = heatmap!(ax3, x, y, diff; colormap=Reverse(:RdBu), colorrange=(-clim_diff, clim_diff))
    Colorbar(fig[row, 6], hm3, width=12)
end

save("comparison_old_vs_new.png", fig, px_per_unit=2)
display(fig)
println("\nSaved → comparison_old_vs_new.png")

end # let

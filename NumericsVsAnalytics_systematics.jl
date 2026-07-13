using CairoMakie, MathTeXEngine, Printf, Statistics
using LinearAlgebra, ExactFieldSolutions, StaticArrays
using JLD2

#Old analytical solution (Duretz)
#Matrix potentials
ϕ_mat(z,   k, γ̇, ε̇, ζ̇, ηm, κm, ri) =  2*ηm*ζ̇/(κm-1)*z - im/2*ηm*γ̇*z      -   (im*γ̇ + 2*ε̇) * ηm * k[1] * ri^2 / z^1
ϕ′_mat(z,  k, γ̇, ε̇, ζ̇, ηm, κm, ri) =  2*ηm*ζ̇/(κm-1)   - im/2*ηm*γ̇        +   (im*γ̇ + 2*ε̇) * ηm * k[1] * ri^2 / z^2
ϕ′′_mat(z, k, γ̇, ε̇, ζ̇, ηm, κm, ri) =                  - 2*(im*γ̇ + 2*ε̇) * ηm * k[1] * ri^2 / z^3
ψ_mat(z,   k, γ̇, ε̇, ζ̇, ηm, κm, ri) =                   (im*γ̇ - 2*ε̇)*ηm*z -   (im*γ̇ + 2*ε̇) * ηm * k[1] * ri^4 / z^3
ψ′_mat(z,  k, γ̇, ε̇, ζ̇, ηm, κm, ri) =                   (im*γ̇ - 2*ε̇)*ηm   + 3*(im*γ̇ + 2*ε̇) * ηm * k[1] * ri^4 / z^4
#Inclusion potentials
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
    ηm, ηi, ξm, ξi, ri, γ̇, ε̇, ζ̇ = params
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

#New analytical solution (Moutzouris_2026)
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

#Parameter sets
function test_params(test)
    ri = 0.1
    if test == :Schmid2003
        ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e10, 1e10
        ε̇ = -1e0;  γ̇ = 0.0;  ζ̇ = 0.0
    elseif test == :Duretz2026_1
        ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e0, 1e0
        ε̇ = -1e0;  γ̇ = 0.0;  ζ̇ = 0.0
    elseif test == :Duretz2026_2
        ηm, ηi = 1.0, 1e-1;  ξm, ξi = 1e0, 1e0
        ε̇ = 1e-10;  γ̇ = 1.0;  ζ̇ = 1e0
    elseif test == :Duretz2026_3
        ηm, ηi = 1.0, 1e-1;  ξm, ξi = 1.0, 1.0
        ε̇ = 1e-10;  γ̇ = 0.0;  ζ̇ = 1e0
    else
        error("Unknown test: $test")
    end
    return (ηm=ηm, ηi=ηi, ξm=ξm, ξi=ξi, ri=ri, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇)
end

#Numerical pseudotransient solver
@views av4_harm(A) = 1.0 ./ (0.25 .* (1.0./A[1:end-1,1:end-1] .+ 1.0./A[2:end,1:end-1] .+ 1.0./A[1:end-1,2:end] .+ 1.0./A[2:end,2:end]))

function Gershgorin_Stokes2D_SchurComplement(ηc, ηv, γ, Δx, Δy, ncx, ncy)
    ηN  = ones(ncx-1, ncy);  ηS  = ones(ncx-1, ncy)
    ηN[:,1:end-1] .= ηv[2:end-1,2:end-1];  ηS[:,2:end] .= ηv[2:end-1,2:end-1]
    ηW  = ηc[1:end-1,:];  ηE = ηc[2:end,:]
    ebW = γ[1:end-1,:];   ebE = γ[2:end,:]
    Cxx = ones(ncx-1, ncy);  Cxy = ones(ncx-1, ncy)
    @. Cxx = abs.(ηN./Δy^2) + abs.(ηS./Δy^2) + abs.(ebE./Δx^2 + (4//3)*ηE./Δx^2) + abs.(ebW./Δx^2 + (4//3)*ηW./Δx^2) + abs.(-(-ηN./Δy - ηS./Δy)./Δy + (ebE./Δx + ebW./Δx)./Δx + ((4//3)*ηE./Δx + (4//3)*ηW./Δx)./Δx)
    @. Cxy = abs.(ebE./(Δx*Δy) - 2//3*ηE./(Δx*Δy) + ηN./(Δx*Δy)) + abs.(ebE./(Δx*Δy) - 2//3*ηE./(Δx*Δy) + ηS./(Δx*Δy)) + abs.(ebW./(Δx*Δy) + ηN./(Δx*Δy) - 2//3*ηW./(Δx*Δy)) + abs.(ebW./(Δx*Δy) + ηS./(Δx*Δy) - 2//3*ηW./(Δx*Δy))
    Dx  = ones(ncx-1, ncy)
    @. Dx .= -(-ηN./Δy - ηS./Δy)./Δy + (ebE./Δx + ebW./Δx)./Δx + ((4//3)*ηE./Δx + (4//3)*ηW./Δx)./Δx

    ηE2 = ones(ncx, ncy-1);  ηW2 = ones(ncx, ncy-1)
    ηE2[1:end-1,:] .= ηv[2:end-1,2:end-1];  ηW2[2:end,:] .= ηv[2:end-1,2:end-1]
    ηS2 = ηc[:,1:end-1];  ηN2 = ηc[:,2:end]
    ebS = γ[:,1:end-1];   ebN = γ[:,2:end]
    Cyy = ones(ncx, ncy-1);  Cyx = ones(ncx, ncy-1)
    @. Cyy = abs.(ηE2./Δx^2) + abs.(ηW2./Δx^2) + abs.(ebN./Δy^2 + (4//3)*ηN2./Δy^2) + abs.(ebS./Δy^2 + (4//3)*ηS2./Δy^2) + abs.((ebN./Δy + ebS./Δy)./Δy + ((4//3)*ηN2./Δy + (4//3)*ηS2./Δy)./Δy - (-ηE2./Δx - ηW2./Δx)./Δx)
    @. Cyx = abs.(ebN./(Δx*Δy) + ηE2./(Δx*Δy) - 2//3*ηN2./(Δx*Δy)) + abs.(ebN./(Δx*Δy) - 2//3*ηN2./(Δx*Δy) + ηW2./(Δx*Δy)) + abs.(ebS./(Δx*Δy) + ηE2./(Δx*Δy) - 2//3*ηS2./(Δx*Δy)) + abs.(ebS./(Δx*Δy) - 2//3*ηS2./(Δx*Δy) + ηW2./(Δx*Δy))
    Dy  = ones(ncx, ncy-1)
    @. Dy .= (ebN./Δy + ebS./Δy)./Δy + ((4//3)*ηN2./Δy + (4//3)*ηS2./Δy)./Δy - (-ηE2./Δx - ηW2./Δx)./Δx

    λmaxVx = 1.0./Dx .* (Cxx .+ Cxy)
    λmaxVy = 1.0./Dy .* (Cyx .+ Cyy)
    return Dx, Dy, λmaxVx, λmaxVy
end
#Stokes solver
@views function Stokes2D(n, test; formulation=:new)

    f_anal = formulation == :new ? Analytics_new : Analytics_old

    Lx, Ly   = 1.0, 1.0
    comp     = true
    one_iter = false

    params = test_params(test)
    ηm, ηi, ξm, ξi, ri, γ̇, ε̇, ζ̇ = params

    # For BCs use the selected analytical solution with physical compressibility
    bc_params = comp ? (ηm=ηm, ηi=ηi, ξm=ξm, ξi=ξm, ri=ri, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇) :
                       (ηm=ηm, ηi=ηi, ξm=1e100*ξm, ξi=1e100*ξm, ri=ri, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇)

    ncx, ncy = n, n
    ϵ        = 1e-7
    iterMax  = 1e5
    nout     = 100
    c_fact   = 0.5
    dτ_local = true
    γfact    = 60
    rel_drop = 1e-3

    Δx, Δy  = Lx/ncx, Ly/ncy
    Pt       = zeros(ncx,   ncy  )
    ∇V       = zeros(ncx,   ncy  )
    Vx       = zeros(ncx+1, ncy+2)
    Vy       = zeros(ncx+2, ncy+1)
    dVx      = zeros(ncx-1, ncy  )
    dVy      = zeros(ncx,   ncy-1)
    Exx      = zeros(ncx,   ncy  )
    Eyy      = zeros(ncx,   ncy  )
    Exy      = zeros(ncx+1, ncy+1)
    Txx      = zeros(ncx,   ncy  )
    Tyy      = zeros(ncx,   ncy  )
    Txy      = zeros(ncx+1, ncy+1)
    Sxx      = zeros(ncx,   ncy  )
    Syy      = zeros(ncx,   ncy  )
    Szz      = zeros(ncx,   ncy  )
    Rx       = zeros(ncx-1, ncy  )
    Ry       = zeros(ncx,   ncy-1)
    Rp       = zeros(ncx,   ncy  )
    Rx0      = zeros(ncx-1, ncy  )
    Ry0      = zeros(ncx,   ncy-1)
    dVxdτ    = zeros(ncx-1, ncy  )
    dVydτ    = zeros(ncx,   ncy-1)
    βVx      = zeros(ncx-1, ncy  )
    βVy      = zeros(ncx,   ncy-1)
    cVx      = zeros(ncx-1, ncy  )
    cVy      = zeros(ncx,   ncy-1)
    αVx      = zeros(ncx-1, ncy  )
    αVy      = zeros(ncx,   ncy-1)
    ηb       = zeros(ncx,   ncy  )
    ηc       = ηm .* ones(ncx,   ncy  )
    ηv       = ηm .* ones(ncx+1, ncy+1)
    ηc_sharp = ηm .* ones(ncx,   ncy  )
    ηv_sharp = ηm .* ones(ncx+1, ncy+1)

    xce, yce = LinRange(-Lx/2-Δx/2, Lx/2+Δx/2, ncx+2), LinRange(-Ly/2-Δy/2, Ly/2+Δy/2, ncy+2)
    xc, yc   = LinRange(-Lx/2+Δx/2, Lx/2-Δx/2, ncx),   LinRange(-Ly/2+Δy/2, Ly/2-Δy/2, ncy)
    xv, yv   = LinRange(-Lx/2, Lx/2, ncx+1),             LinRange(-Ly/2, Ly/2, ncy+1)

    ηv_sharp[(xv).^2 .+ (yv').^2 .< ri^2] .= ηi
    ηc_sharp[(xc).^2 .+ (yc').^2 .< ri^2] .= ηi
    ηc .= av4_harm(ηv_sharp)
    ηv[2:end-1,2:end-1] .= av4_harm(ηc_sharp)
    ηb .= ξm

    γi    = γfact * mean(ηc) .* ones(size(ηc))
    γ_eff = zeros(size(ηb))
    if comp
        γ_num = γi .* ones(size(ηb))
        γ_phy = ηb
        γ_eff = (γ_phy.^-1 .+ γ_num.^-1).^-1
    else
        γ_eff .= γi
    end

    if one_iter
        γ_eff  = ηb
        rel_drop = 1e-11
    end

    Dx, Dy, λmaxVx, λmaxVy = Gershgorin_Stokes2D_SchurComplement(ηc, ηv, γ_eff, Δx, Δy, ncx, ncy)
    if dτ_local
        dτVx = 2.0 ./ sqrt.(λmaxVx) * 0.99
        dτVy = 2.0 ./ sqrt.(λmaxVy) * 0.99
    else
        dτVx = 2.0 / sqrt(maximum(λmaxVx)) * 0.99
        dτVy = 2.0 / sqrt(maximum(λmaxVy)) * 0.99
    end
    βVx .= 2 .* dτVx ./ (2 .+ cVx .* dτVx)
    βVy .= 2 .* dτVy ./ (2 .+ cVy .* dτVy)
    αVx .= (2 .- cVx .* dτVx) ./ (2 .+ cVx .* dτVx)
    αVy .= (2 .- cVy .* dτVy) ./ (2 .+ cVy .* dτVy)

    # Boundary values from analytical solution
    VxS, VxN = zeros(ncx+1), zeros(ncx+1)
    for i in eachindex(VxS)
        VxS[i] = f_anal(@SVector([xv[i]; -Ly/2]); params=bc_params).V[1]
        VxN[i] = f_anal(@SVector([xv[i];  Ly/2]); params=bc_params).V[1]
    end
    VyW, VyE = zeros(ncy+1), zeros(ncy+1)
    for j in eachindex(VyW)
        VyW[j] = f_anal(@SVector([-Lx/2; yv[j]]); params=bc_params).V[2]
        VyE[j] = f_anal(@SVector([ Lx/2; yv[j]]); params=bc_params).V[2]
    end

    # Initial condition: zero interior first (creates non-zero pressure residual),
    # then overwrite all of Vx/Vy from the analytical solution.
    Vx[2:end-1,:] .= 0
    Vy[:,2:end-1] .= 0
    Va = (x=zeros(ncx+1,ncy+2), y=zeros(ncx+2,ncy+1))
    for I in CartesianIndices(Vx)
        i, j      = I[1], I[2]
        sol       = f_anal(@SVector([xv[i]; yce[j]]); params=bc_params)
        Vx[i,j]   = sol.V[1]
        Va.x[i,j] = sol.V[1]
    end
    for I in CartesianIndices(Vy)
        i, j      = I[1], I[2]
        sol       = f_anal(@SVector([xce[i]; yv[j]]); params=bc_params)
        Vy[i,j]   = sol.V[2]
        Va.y[i,j] = sol.V[2]
    end

    # Analytical reference on cell-centre / vertex grids
    Pa = zeros(ncx, ncy)
    Sa = (xx=zeros(ncx,ncy), yy=zeros(ncx,ncy), zz=zeros(ncx,ncy), xy=zeros(ncx+1,ncy+1))
    for I in CartesianIndices(Pa)
        i, j = I[1], I[2]
        sol = f_anal(@SVector([xc[i]; yc[j]]); params=bc_params)
        Pa[i,j]    = sol.p
        Sa.xx[i,j] = sol.τ[1,1] - sol.p
        Sa.yy[i,j] = sol.τ[2,2] - sol.p
        Sa.zz[i,j] = sol.τzz    - sol.p
    end
    for I in CartesianIndices(Sa.xy)
        i, j = I[1], I[2]
        sol = f_anal(@SVector([xv[i]; yv[j]]); params=bc_params)
        Sa.xy[i,j] = sol.τ[1,2]
    end

    # Iteration loop
    errVx0=1.0; errVy0=1.0; errPt0=1.0
    errVx00=1.0; errVy00=1.0
    iter=1; err=2*ϵ
    @time for itPH = 1:50
        Vx[:,1] .= 2*VxS .- Vx[:,2];  Vx[:,end] .= 2*VxN .- Vx[:,end-1]
        Vy[1,:] .= 2*VyW .- Vy[2,:];  Vy[end,:] .= 2*VyE .- Vy[end-1,:]
        ∇V  .= (Vx[2:end,2:end-1] .- Vx[1:end-1,2:end-1])./Δx .+ (Vy[2:end-1,2:end] .- Vy[2:end-1,1:end-1])./Δy
        Exx .= (Vx[2:end,2:end-1] .- Vx[1:end-1,2:end-1])./Δx .- 1/3 .* ∇V
        Eyy .= (Vy[2:end-1,2:end] .- Vy[2:end-1,1:end-1])./Δy .- 1/3 .* ∇V
        Exy .= 0.5 .* ((Vx[:,2:end] .- Vx[:,1:end-1])./Δy .+ (Vy[2:end,:] .- Vy[1:end-1,:])./Δx)
        Txx .= 2 .* ηc .* Exx
        Tyy .= 2 .* ηc .* Eyy
        Txy .= 2 .* ηv .* Exy
        Rx  .= (.-(Pt[2:end,:] .- Pt[1:end-1,:])./Δx .+ (Txx[2:end,:] .- Txx[1:end-1,:])./Δx .+ (Txy[2:end-1,2:end] .- Txy[2:end-1,1:end-1])./Δy)
        Ry  .= (.-(Pt[:,2:end] .- Pt[:,1:end-1])./Δy .+ (Tyy[:,2:end] .- Tyy[:,1:end-1])./Δy .+ (Txy[2:end,2:end-1] .- Txy[1:end-1,2:end-1])./Δx)
        Rp  .= .-∇V .- comp*Pt./ηb
        errVx = norm(Rx); errVy = norm(Ry); errPt = norm(Rp)
        if itPH==1 errVx0=errVx; errVy0=errVy; errPt0=errPt end
        err = maximum([errVx/errVx0, errVy/errVy0, errPt/errPt0])
        @printf("itPH=%02d  err=%1.3e  [Rx=%1.3e  Ry=%1.3e  Rp=%1.3e]\n", itPH, err, errVx/errVx0, errVy/errVy0, errPt/errPt0)
        if err < ϵ break end
        ϵ_vel = err * rel_drop
        itPT  = 0
        while err > ϵ_vel && itPT <= iterMax
            itPT += 1; itg = iter
            Rx0  .= Rx;  Ry0 .= Ry
            Vx[:,1] .= 2*VxS .- Vx[:,2];  Vx[:,end] .= 2*VxN .- Vx[:,end-1]
            Vy[1,:] .= 2*VyW .- Vy[2,:];  Vy[end,:] .= 2*VyE .- Vy[end-1,:]
            ∇V  .= (Vx[2:end,2:end-1] .- Vx[1:end-1,2:end-1])./Δx .+ (Vy[2:end-1,2:end] .- Vy[2:end-1,1:end-1])./Δy
            Exx .= (Vx[2:end,2:end-1] .- Vx[1:end-1,2:end-1])./Δx .- 1/3 .* ∇V
            Eyy .= (Vy[2:end-1,2:end] .- Vy[2:end-1,1:end-1])./Δy .- 1/3 .* ∇V
            Exy .= 0.5 .* ((Vx[:,2:end] .- Vx[:,1:end-1])./Δy .+ (Vy[2:end,:] .- Vy[1:end-1,:])./Δx)
            Rp  .= .-∇V .- comp*Pt./ηb
            Txx .= 2 .* ηc .* Exx .- γ_eff .* Rp
            Tyy .= 2 .* ηc .* Eyy .- γ_eff .* Rp
            Txy .= 2 .* ηv .* Exy
            Rx  .= (1 ./ Dx) .* (.-(Pt[2:end,:] .- Pt[1:end-1,:])./Δx .+ (Txx[2:end,:] .- Txx[1:end-1,:])./Δx .+ (Txy[2:end-1,2:end] .- Txy[2:end-1,1:end-1])./Δy)
            Ry  .= (1 ./ Dy) .* (.-(Pt[:,2:end] .- Pt[:,1:end-1])./Δy .+ (Tyy[:,2:end] .- Tyy[:,1:end-1])./Δy .+ (Txy[2:end,2:end-1] .- Txy[1:end-1,2:end-1])./Δx)
            dVxdτ .= αVx .* dVxdτ .+ Rx
            dVydτ .= αVy .* dVydτ .+ Ry
            Vx[2:end-1,2:end-1] .+= dVxdτ .* βVx .* dτVx
            Vy[2:end-1,2:end-1] .+= dVydτ .* βVy .* dτVy
            if mod(iter, nout) == 0
                errVx = norm(Dx .* Rx); errVy = norm(Dy .* Ry)
                if iter == nout errVx00=errVx; errVy00=errVy end
                err = maximum([min(errVx, errVx/errVx00), min(errVy, errVy/errVy00)])
                push!([], errVx/errVx00)
                dVx .= dVxdτ .* βVx .* dτVx
                dVy .= dVydτ .* βVy .* dτVy
                λminV  = abs(sum(dVx .* Dx .* (Rx .- Rx0)) + sum(dVy .* Dy .* (Ry .- Ry0))) / (sum(dVx .* Dx .* dVx) + sum(dVy .* Dy .* dVy))
                cVx .= 2 * sqrt.(λminV) * c_fact * 0.5
                cVy .= 2 * sqrt.(λminV) * c_fact * 0.5
                βVx .= 2 .* dτVx ./ (2 .+ cVx .* dτVx)
                βVy .= 2 .* dτVy ./ (2 .+ cVy .* dτVy)
                αVx .= (2 .- cVx .* dτVx) ./ (2 .+ cVx .* dτVx)
                αVy .= (2 .- cVy .* dτVy) ./ (2 .+ cVy .* dτVy)
            end
            iter += 1
        end
        Pt .+= γ_eff .* Rp
    end

    Sxx = -Pt .+ Txx
    Syy = -Pt .+ Tyy
    Szz = -Pt .+ (-Txx .- Tyy)

    V = (x=Vx, y=Vy)
    S = (xx=Sxx, yy=Syy, zz=Szz, xy=Txy)

    # Errors
    Ve = ( 
        x  = abs.(V.x .- Va.x),
        y  = abs.(V.y .- Va.y),
    )
    Pe  = abs.(Pt .- Pa)
    Se  = (
        xx = abs.(S.xx .- Sa.xx),
        yy = abs.(S.yy .- Sa.yy),
        zz = abs.(S.zz .- Sa.zz),
        xy = abs.(S.xy .- Sa.xy),
    )
    @printf("mean|dVx| = %1.3e   mean|dVy| = %1.3e   mean|dP| = %1.3e   
        mean|dSxx| = %1.3e   mean|dSyy| = %1.3e   mean|dSzz| = %1.3e   mean|dSzz| = %1.3e\n",
        mean(abs, Ve.x), mean(abs, Ve.y), mean(abs, Pe),
        mean(abs, Se.xx), mean(abs, Se.yy), mean(abs, Se.zz), mean(abs, Se.xy))
    
    L1 = (Vx=mean(Ve.x), Vy=mean(Ve.y), P=mean(Pe), σxx=mean(Se.xx), σyy=mean(Se.yy), σzz=mean(Se.zz), σxy=mean(Se.xy))

    @show iter/ncx
    return L1, V, Pt, S, Va, Pa, Sa, xv, yv, xce, yce, xc, yc
end

# Compare solutions
let
    formulation = :new             #:old | :new
    Lx, Ly = 1.0, 1.0
    nc          = [51, 101, 201, 401]              #resolution

    tests      = (:Schmid2003, :Duretz2026_1, :Duretz2026_2, :Duretz2026_3)

    L1  = (
        h   = zeros(length(tests), length(nc)), 
        Vx  = zeros(length(tests), length(nc)), 
        Vy  = zeros(length(tests), length(nc)), 
        P   = zeros(length(tests), length(nc)), 
        σxx = zeros(length(tests), length(nc)), 
        σyy = zeros(length(tests), length(nc)), 
        σzz = zeros(length(tests), length(nc)), 
        σxy = zeros(length(tests), length(nc))
    )

    for itest in eachindex(tests), ires in eachindex(nc)

        errors, V, P, S, Va, Pa, Sa, xv, yv, xce, yce, xc, yc =
        Stokes2D(nc[ires], tests[itest]; formulation=formulation)
     
        L1.h[itest, ires]   = Lx/nc[ires]
        L1.Vx[itest, ires]  = errors.Vx  
        L1.Vy[itest, ires]  = errors.Vy  
        L1.P[itest, ires]   = errors.P   
        L1.σxx[itest, ires] = errors.σxx 
        L1.σyy[itest, ires] = errors.σyy 
        L1.σzz[itest, ires] = errors.σzz 
        L1.σxy[itest, ires] = errors.σxy
    end

    jldsave("TruncationSystematics.jld2", errors=L1)
end

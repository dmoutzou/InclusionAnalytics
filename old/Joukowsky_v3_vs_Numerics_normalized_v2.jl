#Dimitrios Moutzouris, 2026, Analytical solution for compressible elliptical inclusions in host #
#This is a comparison file between analytics and numerics
#
# Normalised version: the numerical domain is the fixed unit box [-1,1]^2.
# Physical coordinate  X_phys = S * x_norm,   with  S = scale * a.
# The analytical solution lives in physical space, so we evaluate it at S*x_norm
# and rescale the fields:  V_num = V_phys / S ;  P and tau are scale-invariant.
using CairoMakie, MathTeXEngine, Printf, Statistics, LinearAlgebra, ExactFieldSolutions, StaticArrays

@views av4_harm(A) = 1.0./( 0.25.*(1.0./A[1:end-1,1:end-1].+1.0./A[2:end,1:end-1].+1.0./A[1:end-1,2:end].+1.0./A[2:end,2:end]) )
include("Joukowsky_v4_DM_function.jl") 

# Analytical field at PHYSICAL coordinate
f_anal(X; params) = Analytics_new(X; params=params)

function preset(test)
    if test == :Schmid2003
        ηm = 1.0; ηi = 1000.0;  ξm, ξi = 1e10, 1e10
        ε̇ = -0.5;  γ̇ = 0.0;  ζ̇ = 0.0;  t = 2.0;  α = 0.0; ri = 0.2
    elseif test == :Duretz2026_1
        ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e0, 1e0
        ε̇ =-1e0;  γ̇ = 0.0;  ζ̇ = 0.0;  t = 2.0;  α = 0.0; ri = 0.2
    elseif test == :Duretz2026_2
        ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e0, 10e0
        ε̇ = 1e-10;  γ̇ = 0.0;  ζ̇ = 1e0;  t = 2.0;  α = 0.0; ri = 0.2
    elseif test == :Duretz2026_3
        ηm, ηi = 1.0, 1000.0;  ξm, ξi = 1e10, 1e10
        ε̇ = 0.0;  γ̇ = 1.0;  ζ̇ = 0.0;  t = 2.0;  α = 30/180*pi; ri = 0.2
    else
        error("unknown test $test")
    end
    a, b     = ellipse_axes(t)   # true physical semi-axes (a >= 2 always)
    sc       = a / ri 
    return (ηm=ηm, ηi=ηi, ξm=ξm, ξi=ξi, t=t, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇, α=α, ri=ri, sc=sc)
end

# Can be replaced by AD
function Gershgorin_Stokes2D_SchurComplement(ηc, ηv, γ, Δx, Δy, ncx  ,ncy)

    ηN    = ones(ncx-1, ncy)
    ηS    = ones(ncx-1, ncy)
    ηN[:,1:end-1] .= ηv[2:end-1,2:end-1]
    ηS[:,2:end-0] .= ηv[2:end-1,2:end-1]
    ηW    = ηc[1:end-1,:]
    ηE    = ηc[2:end-0,:]
    ebW   = γ[1:end-1,:]
    ebE   = γ[2:end-0,:]
    Cxx   = ones(ncx-1, ncy)
    Cxy   = ones(ncx-1, ncy)
    @. Cxx = abs.(ηN ./ Δy .^ 2) + abs.(ηS ./ Δy .^ 2) + abs.(ebE ./ Δx .^ 2 + (4 // 3) * ηE ./ Δx .^ 2) + abs.(ebW ./ Δx .^ 2 + (4 // 3) * ηW ./ Δx .^ 2) + abs.(-(-ηN ./ Δy - ηS ./ Δy) ./ Δy + (ebE ./ Δx + ebW ./ Δx) ./ Δx + ((4 // 3) * ηE ./ Δx + (4 // 3) * ηW ./ Δx) ./ Δx)
    @. Cxy = abs.(ebE ./ (Δx .* Δy) - 2 // 3 * ηE ./ (Δx .* Δy) + ηN ./ (Δx .* Δy)) + abs.(ebE ./ (Δx .* Δy) - 2 // 3 * ηE ./ (Δx .* Δy) + ηS ./ (Δx .* Δy)) + abs.(ebW ./ (Δx .* Δy) + ηN ./ (Δx .* Δy) - 2 // 3 * ηW ./ (Δx .* Δy)) + abs.(ebW ./ (Δx .* Δy) + ηS ./ (Δx .* Δy) - 2 // 3 * ηW ./ (Δx .* Δy))

    Dx  = ones(ncx-1, ncy)
    @. Dx .= -(-ηN ./ Δy - ηS ./ Δy) ./ Δy + (ebE ./ Δx + ebW ./ Δx) ./ Δx + ((4 // 3) * ηE ./ Δx + (4 // 3) * ηW ./ Δx) ./ Δx

    ηE    = ones(ncx, ncy-1)
    ηW    = ones(ncx, ncy-1)
    ηE[1:end-1,:] .= ηv[2:end-1,2:end-1]
    ηW[2:end-0,:] .= ηv[2:end-1,2:end-1]
    ηS    = ηc[:,1:end-1]
    ηN    = ηc[:,2:end-0]
    ebS  = γ[:,1:end-1]
    ebN  = γ[:,2:end-0]
    Cyy  = ones(ncx, ncy-1)
    Cyx  = ones(ncx, ncy-1)
    @. Cyy = abs.(ηE ./ Δx .^ 2) + abs.(ηW ./ Δx .^ 2) + abs.(ebN ./ Δy .^ 2 + (4 // 3) * ηN ./ Δy .^ 2) + abs.(ebS ./ Δy .^ 2 + (4 // 3) * ηS ./ Δy .^ 2) + abs.((ebN ./ Δy + ebS ./ Δy) ./ Δy + ((4 // 3) * ηN ./ Δy + (4 // 3) * ηS ./ Δy) ./ Δy - (-ηE ./ Δx - ηW ./ Δx) ./ Δx)
    @. Cyx = abs.(ebN ./ (Δx .* Δy) + ηE ./ (Δx .* Δy) - 2 // 3 * ηN ./ (Δx .* Δy)) + abs.(ebN ./ (Δx .* Δy) - 2 // 3 * ηN ./ (Δx .* Δy) + ηW ./ (Δx .* Δy)) + abs.(ebS ./ (Δx .* Δy) + ηE ./ (Δx .* Δy) - 2 // 3 * ηS ./ (Δx .* Δy)) + abs.(ebS ./ (Δx .* Δy) - 2 // 3 * ηS ./ (Δx .* Δy) + ηW ./ (Δx .* Δy))

    Dy  = ones(ncx, ncy-1)
    @. Dy .= (ebN ./ Δy + ebS ./ Δy) ./ Δy + ((4 // 3) * ηN ./ Δy + (4 // 3) * ηS ./ Δy) ./ Δy - (-ηE ./ Δx - ηW ./ Δx) ./ Δx

    λmaxVx = 1.0./Dx .* (Cxx .+ Cxy)
    λmaxVy = 1.0./Dy .* (Cyx .+ Cyy)

    return Dx, Dy, λmaxVx, λmaxVy
end


#2D Stokes solve on the NORMALISED unit box [-1,1]^2.
# `scale` is the old `zoom` knob: physical half-width = scale * a, so the
# normalised ellipse has semi-axes (a/S, b/S) = (1/scale, b/(scale*a)).
# scale > 1 keeps the inclusion inside the box (scale = 1 => tips touch the walls).
@views function Stokes2D(n, test; scale = 3.0)

    #Pull the same physical parameters
    params   = preset(test)
    ηm       = params.ηm
    γ̇, ε̇, ζ̇  = params.γ̇, params.ε̇, params.ζ̇
    a, b     = ellipse_axes(params.t)   # true physical semi-axes (a >= 2 always)
    comp     = true         # keep true :D
    one_iter = false        # true: simple DR --- false: PH/DR

    # (Pseudo-)compressibility: optionally push toward the incompressible limit
    if !comp
        params = merge(params, (ξm = 1e100*params.ξm, ξi = 1e100*params.ξi))
    end

    # --- normalisation --------------------------------------------------
    # X_phys = S * x_norm,  x_norm ∈ [-1,1]
    S        = a / params.ri              # physical half-width of the box
    an, bn   = a/S, b/S              # normalised semi-axes (an = 1/scale)
    Lx, Ly   = 1.0, 1.0             # fixed unit box [-1,1]^2
    # analytics evaluated at a NORMALISED point; returns PHYSICAL-unit fields.
    # Remember to divide velocities by S; p and tau are scale-invariant.
    anal(x, y) = f_anal(@SVector([x, y]); params=params)
    # --------------------------------------------------------------------

    # Numerics
    ncx, ncy = n, n         # numerical grid resolution
    ϵ        = 1e-7         # tolerance
    iterMax  = 1e5          # max number of iters
    nout     = 100          # check frequency
    c_fact   = 0.5          # damping factor
    C        = 0.99
    dτ_local = false        # helps a little bit sometimes, sometimes not!
    γfact    = 5            # penalty: multiplier to the arithmetic mean of η
    rel_drop = 1e-3         # relative drop of velocity residual per PH iteration
    nPH      = 50
    # Preprocessing
    Δx, Δy  = Lx/ncx, Ly/ncy
    # Array initialisation
    Pt       = zeros(ncx  ,ncy  )
    ∇V       = zeros(ncx  ,ncy  )
    Vx       = zeros(ncx+1,ncy+2)
    Vy       = zeros(ncx+2,ncy+1)
    dVx      = zeros(ncx-1,ncy  )
    dVy      = zeros(ncx  ,ncy-1)
    Exx      = zeros(ncx  ,ncy  )
    Eyy      = zeros(ncx  ,ncy  )
    Exy      = zeros(ncx+1,ncy+1)
    Txx      = zeros(ncx  ,ncy  )
    Tyy      = zeros(ncx  ,ncy  )
    Txy      = zeros(ncx+1,ncy+1)
    Sxx      = zeros(ncx  ,ncy  )
    Syy      = zeros(ncx  ,ncy  )
    Szz      = zeros(ncx  ,ncy  )
    Rx       = zeros(ncx-1,ncy  )
    Ry       = zeros(ncx  ,ncy-1)
    Rp       = zeros(ncx  ,ncy  )
    Rx0      = zeros(ncx-1,ncy  )
    Ry0      = zeros(ncx  ,ncy-1)
    dVxdτ    = zeros(ncx-1,ncy  )
    dVydτ    = zeros(ncx  ,ncy-1)
    βVx      = zeros(ncx-1,ncy  )
    βVy      = zeros(ncx  ,ncy-1)
    cVx      = zeros(ncx-1,ncy  )  # this disappears is dτ is not local
    cVy      = zeros(ncx  ,ncy-1)  # this disappears is dτ is not local
    αVx      = zeros(ncx-1,ncy  )  # this disappears is dτ is not local
    αVy      = zeros(ncx  ,ncy-1)  # this disappears is dτ is not local
    ηb       = zeros(ncx  ,ncy  )
    ηc       = ηm .* ones(ncx  ,ncy  )
    ηv       = ηm .* ones(ncx+1,ncy+1)
    ηc_sharp = ηm .* ones(ncx  ,ncy  )
    ηv_sharp = ηm .* ones(ncx+1,ncy+1)
    # Initialisation (normalised grid on [-1,1])
    xce, yce = LinRange(-Lx/2-Δx/2, Lx/2+Δx/2, ncx+2), LinRange(-Ly/2-Δy/2, Ly/2+Δy/2, ncy+2)
    xc, yc   = LinRange(-Lx/2+Δx/2, Lx/2-Δx/2, ncx), LinRange(-Ly/2+Δy/2, Ly/2-Δy/2, ncy)
    xv, yv   = LinRange(-Lx/2, Lx/2, ncx+1), LinRange(-Ly/2, Ly/2, ncy+1)
    ηv_sharp   .= ηm
    ηv_sharp[(xv./an).^2 .+ (yv'./bn).^2 .< 1.0 ] .= params.ηi
    ηc_sharp   .= ηm
    ηc_sharp[(xc./an).^2 .+ (yc'./bn).^2 .< 1.0 ] .= params.ηi
    # Use sharp viscosity
    ηc .= ηc_sharp
    ηv .= ηv_sharp
    # Harmonic averaging mimicking PIC interpolation
    ηc    .= av4_harm(ηv_sharp)
    ηv[2:end-1,2:end-1] .= av4_harm(ηc_sharp)
    # Set viscosity
    ηb    .= params.ξm
    ηb[(xc./an).^2 .+ (yc'./bn).^2 .< 1.0 ] .= params.ξi

    # Select γ
    γi   = γfact*mean(ηc).*ones(size(ηc))

    # (Pseudo-)compressibility
    γ_eff = zeros(size(ηb))
    if comp
        γ_num = γi.*ones(size(ηb)) * 1.0
        γ_phy = ηb
        γ_eff = (γ_phy.^-1 .+ γ_num.^-1).^-1
    else
        γ_eff .= γi
        γ_eff .= γ_eff
    end

    ##########################################
    if one_iter
        # Converges in on PH iter
        γ_phy = ηb
        γ_eff = ηb
        γ_num = ηb
        rel_drop = 1e-11
    end
    ##########################################

    # Optimal pseudo-time steps - can be replaced by AD
    Dx, Dy, λmaxVx, λmaxVy = Gershgorin_Stokes2D_SchurComplement(ηc, ηv, γ_eff, Δx, Δy, ncx ,ncy)
    # Select dτ
    if dτ_local
        dτVx =  2.0./sqrt.(λmaxVx)*C
        dτVy =  2.0./sqrt.(λmaxVy)*C
    else
        dτVx =  2.0./sqrt.(maximum(λmaxVx))*C
        dτVy =  2.0./sqrt.(maximum(λmaxVy))*C
    end
    βVx .= 2 .* dτVx ./ (2 .+ cVx.*dτVx)
    βVy .= 2 .* dτVy ./ (2 .+ cVy.*dτVy)
    αVx .= (2 .- cVx.*dτVx) ./ (2 .+ cVx.*dτVx)
    αVy .= (2 .- cVy.*dτVy) ./ (2 .+ cVy.*dτVy)
    # Initial condition (normalised: V_num = ε̇ * x_norm already, no /S needed)
    Vx     .=   ε̇.*xv .+  0*yce' .+  ζ̇ .*xv
    Vy     .=   0*xce .-  ε̇.*yv' .+  ζ̇ .*yv'
    Vx_ini  = copy(Vx)
    Vy_ini  = copy(Vy)
    Va = (
        x = zeros(ncx+1, ncy+2),
        y = zeros(ncx+2, ncy+1)
    )
    @time for I in CartesianIndices(Vx)
        i, j      = I[1], I[2]
        sol       = anal(xv[i], yce[j])
        Vx[i,j]   = sol.V[1]
        Va.x[i,j] = sol.V[1]
    end

    for I in CartesianIndices(Vy)
        i, j      = I[1], I[2]
        sol       = anal(xce[i], yv[j])
        Vy[i,j]   = sol.V[2]
        Va.y[i,j] = sol.V[2]
    end
    Vx[2:end-1,:] .= 0 # ensure non zero initial pressure residual
    Vy[:,2:end-1] .= 0 # ensure non zero initial pressure residual

    #boundary conditions (physical coord in via `anal`, velocity  out)
    VxS, VxN = zeros(size(Vx,1)), zeros(size(Vx,1))
    for i in eachindex(VxS)
        VxS[i] = anal(xv[i], -Ly/2).V[1]
        VxN[i] = anal(xv[i],  Ly/2).V[1]
    end

    VyW, VyE = zeros(size(Vy,2)), zeros(size(Vy,2))
    for j in eachindex(VyW)
        VyW[j] = anal(-Lx/2, yv[j]).V[2]
        VyE[j] = anal( Lx/2, yv[j]).V[2]
    end

    Pa = zeros(ncx, ncy)
    Sa = (
        xx = zeros(ncx, ncy),
        yy = zeros(ncx, ncy),
        zz = zeros(ncx, ncy),
        xy = zeros(ncx+1, ncy+1),
    )
    for I in CartesianIndices(Pa)
        i, j       = I[1], I[2]
        sol        = anal(xc[i], yc[j])   # p, tau are scale-invariant: no /S
        Pa[i,j]    = sol.p
        Sa.xx[i,j] = sol.τ[1,1] - sol.p
        Sa.yy[i,j] = sol.τ[2,2] - sol.p
        Sa.zz[i,j] = sol.τzz    - sol.p
    end

    for I in CartesianIndices(Sa.xy)
        i, j        = I[1], I[2]
        Sa.xy[i,j]  = anal(xv[i], yv[j]).τ[1,2]
    end

    # Iteration loop
    errVx0 = 1.0;  errVy0 = 1.0;  errPt0 = 1.0
    errVx00= 1.0;  errVy00= 1.0;
    iter=1; err=2*ϵ; err_evo_V=[]; err_evo_P=[]; err_evo_it=[]
    @time for itPH = 1:nPH
        # Boundaries
        Vx[:,1] .= 2*VxS .- Vx[:,2]; Vx[:,end] .= 2*VxN .- Vx[:,end-1]
        Vy[1,:] .= 2*VyW .- Vy[2,:]; Vy[end,:] .= 2*VyE .- Vy[end-1,:]
        # Divergence
        ∇V    .= (Vx[2:end,2:end-1] .- Vx[1:end-1,2:end-1])./Δx .+ (Vy[2:end-1,2:end] .- Vy[2:end-1,1:end-1])./Δy
        # Deviatoric strain rate
        Exx   .= (Vx[2:end,2:end-1] .- Vx[1:end-1,2:end-1])./Δx .- 1.0/3.0.*∇V
        Eyy   .= (Vy[2:end-1,2:end] .- Vy[2:end-1,1:end-1])./Δy .- 1.0/3.0.*∇V
        Exy   .= 0.5.*((Vx[:,2:end] .- Vx[:,1:end-1])./Δy .+ (Vy[2:end,:] .- Vy[1:end-1,:])./Δx)
        # Deviatoric stress
        Txx   .= 2.0.*ηc.*Exx
        Tyy   .= 2.0.*ηc.*Eyy
        Txy   .= 2.0.*ηv.*Exy
        # Residuals
        Rx    .= (.-(Pt[2:end,:] .- Pt[1:end-1,:])./Δx .+ (Txx[2:end,:] .- Txx[1:end-1,:])./Δx .+ (Txy[2:end-1,2:end] .- Txy[2:end-1,1:end-1])./Δy)
        Ry    .= (.-(Pt[:,2:end] .- Pt[:,1:end-1])./Δy .+ (Tyy[:,2:end] .- Tyy[:,1:end-1])./Δy .+ (Txy[2:end,2:end-1] .- Txy[1:end-1,2:end-1])./Δx)
        Rp    .= .-∇V .- comp*Pt./ηb
        
        # Residual check
        errVx = norm(Rx); errVy = norm(Ry); errPt = norm(Rp)
        if itPH==1 errVx0=errVx; errVy0=errVy; errPt0=errPt; end
        err = maximum([errVx/errVx0, errVy/errVy0, errPt/errPt0])
        @printf("itPH = %02d iter = %06d iter/nx = %03d, err = %1.3e norm[Rx=%1.3e, Ry=%1.3e, Rp=%1.3e] \n", itPH, iter, iter/ncx, err, errVx/errVx0, errVy/errVy0, errPt/errPt0)
        if (err<ϵ) break end
        # Set tolerance of velocity solve proportional to residual
        ϵ_vel = err*rel_drop
        itPT = 0.
        while (err>ϵ_vel && itPT<=iterMax)
            itPT   += 1
            itg    = iter
            # Pseudo-old dudes
            Rx0   .= Rx
            Ry0   .= Ry
            # Boundaries
            Vx[:,1] .= 2*VxS .- Vx[:,2]; Vx[:,end] .= 2*VxN .- Vx[:,end-1]
            Vy[1,:] .= 2*VyW .- Vy[2,:]; Vy[end,:] .= 2*VyE .- Vy[end-1,:]
            # Vx[:,1] .= Vx[:,2]; Vx[:,end] .= Vx[:,end-1]
            # Vy[1,:] .= Vy[2,:]; Vy[end,:] .= Vy[end-1,:]
            # Divergence
            ∇V    .= (Vx[2:end,2:end-1] .- Vx[1:end-1,2:end-1])./Δx .+ (Vy[2:end-1,2:end] .- Vy[2:end-1,1:end-1])./Δy
            # Deviatoric strain rate
            Exx   .= (Vx[2:end,2:end-1] .- Vx[1:end-1,2:end-1])./Δx .- 1.0/3.0.*∇V
            Eyy   .= (Vy[2:end-1,2:end] .- Vy[2:end-1,1:end-1])./Δy .- 1.0/3.0.*∇V
            Exy   .= 0.5.*((Vx[:,2:end] .- Vx[:,1:end-1])./Δy .+ (Vy[2:end,:] .- Vy[1:end-1,:])./Δx)
            # "Deviatoric" stress
            Rp    .= .-∇V .- comp*Pt./ηb
            Txx   .= 2.0.*ηc.*Exx .- γ_eff .* Rp
            Tyy   .= 2.0.*ηc.*Eyy .- γ_eff .* Rp
            Txy   .= 2.0.*ηv.*Exy
            # Residuals
            Rx    .= (1.0./Dx).*(.-(Pt[2:end,:] .- Pt[1:end-1,:])./Δx .+ (Txx[2:end,:] .- Txx[1:end-1,:])./Δx .+ (Txy[2:end-1,2:end] .- Txy[2:end-1,1:end-1])./Δy)
            Ry    .= (1.0./Dy).*(.-(Pt[:,2:end] .- Pt[:,1:end-1])./Δy .+ (Tyy[:,2:end] .- Tyy[:,1:end-1])./Δy .+ (Txy[2:end,2:end-1] .- Txy[1:end-1,2:end-1])./Δx)
            # Damping-pong
            dVxdτ .= αVx.*dVxdτ .+ Rx
            dVydτ .= αVy.*dVydτ .+ Ry
            # PT updates
            Vx[2:end-1,2:end-1] .+= dVxdτ.*βVx.*dτVx
            Vy[2:end-1,2:end-1] .+= dVydτ.*βVy.*dτVy
            # Residual check
            if mod(iter, nout)==0
                errVx = norm(Dx.*Rx); errVy = norm(Dy.*Ry)
                if iter==nout errVx00=errVx; errVy00=errVy; end
                err = maximum([min(errVx, errVx./errVx00), min(errVy, errVy./errVy00)])
                push!(err_evo_V, errVx/errVx00); push!(err_evo_P, errPt/errPt0); push!(err_evo_it, itg)
                dVx .= dVxdτ.*βVx.*dτVx
                dVy .= dVydτ.*βVy.*dτVy
                # @printf("iter = %d, err = %1.3e norm[Rx=%1.3e, Ry=%1.3e] \n", iter, err, errVx, errVy)
                # λminV  = abs.((sum(dVx.*(Rx .- Rx0))) + abs.((sum(dVy.*(Ry .- Ry0))) )/ ( sum(dVx.*dVx)) + sum(dVy.*dVy) )
                λminV  = abs(  sum(dVx.*Dx.*(Rx .- Rx0)) + sum(dVy.*Dy.*(Ry .- Ry0))  ) / (sum(dVx.*Dx.*dVx) .+ sum(dVy.*Dy.*dVy))
                cVx .= 2*sqrt.(λminV)*c_fact*0.5
                cVy .= 2*sqrt.(λminV)*c_fact*0.5
                βVx .= 2 .* dτVx ./ (2 .+ cVx.*dτVx)
                βVy .= 2 .* dτVy ./ (2 .+ cVy.*dτVy)
                αVx .= (2 .- cVx.*dτVx) ./ (2 .+ cVx.*dτVx)
                αVy .= (2 .- cVy.*dτVy) ./ (2 .+ cVy.*dτVy)
            end
            iter += 1
        end
        Pt .+= γ_eff.*Rp
    end
    Sxx .= .-Pt .+ Txx
    Syy .= .-Pt .+ Tyy
    Szz .= .-Pt .+ (.-Txx .- Tyy)
    Ve = (
        x  = abs.(Vx .- Va.x),
        y  = abs.(Vy .- Va.y),
    )
    Pe  = abs.(Pt .- Pa)
    Se  = (
        xx = abs.(Sxx .- Sa.xx),
        yy = abs.(Syy .- Sa.yy),
        zz = abs.(Szz .- Sa.zz),
        xy = abs.(Txy .- Sa.xy),
    )
    @show iter/ncx
    return (; Vx, Vy, Pt, Va, Pa, Ve, Pe, Se, Lx, Ly, xv, yv, xc, yc, xce, yce, a, b, an, bn, S)
end

let
    test = :Schmid2003
    #test = :Duretz2026_1
    #test = :Duretz2026_2
    #test = :Duretz2026_3

    nc_list = [101]
    scale   = 3.0            # physical half-width = scale * ellipse semi-major axis a

    # Storage for convergence
    err_Vx = Float64[]
    err_Vy = Float64[]
    err_P  = Float64[]

    for nc in nc_list

        num    = Stokes2D(nc, test; scale=scale)
        params = preset(test)

        push!(err_Vx, maximum(abs, num.Ve.x))
        push!(err_Vy, maximum(abs, num.Ve.y))
        push!(err_P,  maximum(abs, num.Pe  ))

        field_data = [
            ("Vx", num.xv,  num.yce, num.Va.x, num.Vx, num.Ve.x,),
            ("Vy", num.xce, num.yv,  num.Va.y, num.Vy, num.Ve.y,),
            ("P",  num.xc,  num.yc,  num.Pa,   num.Pt, num.Pe,  ),
        ]

        fig = Figure(size=(1100, 900), fontsize=12)
        Label(fig[0, 1:2], "Analytical",        tellwidth=false, fontsize=14, font=:bold)
        Label(fig[0, 3:4], "Numerical",         tellwidth=false, fontsize=14, font=:bold)
        Label(fig[0, 5:6], "Difference",        tellwidth=false, fontsize=14, font=:bold)
        Label(fig[-1, 1:6], "nc = $nc  |  test = $test  |  t = $(params.t)", tellwidth=false, fontsize=15, font=:bold)

        for (row, (name, gx, gy, an, numv, diff)) in enumerate(field_data)
            an_c  = an   .- mean(an)
            num_c = numv .- mean(numv)

            ax1 = Axis(fig[row, 1], aspect=DataAspect(), ylabel=name)
            hm1 = heatmap!(ax1, gx, gy, an_c;  colormap=(Reverse(:matter), 1))
            Colorbar(fig[row, 2], hm1, width=12)

            ax2 = Axis(fig[row, 3], aspect=DataAspect())
            hm2 = heatmap!(ax2, gx, gy, num_c; colormap=(Reverse(:matter), 1))
            Colorbar(fig[row, 4], hm2, width=12)

            ax3 = Axis(fig[row, 5], aspect=DataAspect())
            hm3 = heatmap!(ax3, gx, gy, diff;  colormap=(Reverse(:matter), 1))
            Colorbar(fig[row, 6], hm3, width=12)
        end

        save("./figures/comparison_JoukowskyDMvsNumerics_nc$(nc)_$(test).png", fig, px_per_unit=2)
        display(fig)
        println("nc = $nc done  |  max|dVx|=$(err_Vx[end])  max|dVy|=$(err_Vy[end])  max|dP|=$(err_P[end])")
    end

end

#=
# Convergence study: mean error vs resolution for all four tests

let
    tests   = [ :Schmid2003, :Duretz2026_2]
    titles  = ["Test 1", "Test 2"]
    nc_list = [51, 101]
    scale   = 3.0

    fields = [
        (:vx,  "vx",  :steelblue),
        (:vy,  "vy",  :orange),
        (:p,   "p",   :seagreen),
    ]

    fig = Figure(size = (850, 650), fontsize = 12)
    plots_for_legend  = []
    labels_for_legend = String[]

    for (k, test) in enumerate(tests)
        row, col = divrem(k-1, 2) .+ (1, 1)
        ax = Axis(fig[row, col]; xscale = log10, yscale = log10,
                   title = titles[k], xlabel = "1/h", ylabel = "ε")

        invh   = Float64[]
        ε      = Dict(sym => Float64[] for (sym, _, _) in fields)
        params = preset(test)

        for nc in nc_list
            num = Stokes2D(nc, test; scale = scale)
            push!(invh, nc/num.Lx)          # normalised 1/h  (Lx = 2)
            S = num.S

            Vx_an = [f_anal(@SVector([S*num.xv[i],  S*num.yce[j]]); params = params).V[1]/S for i in axes(num.Vx, 1), j in axes(num.Vx, 2)]
            Vy_an = [f_anal(@SVector([S*num.xce[i], S*num.yv[j] ]); params = params).V[2]/S for i in axes(num.Vy, 1), j in axes(num.Vy, 2)]
            P_an  = [f_anal(@SVector([S*num.xc[i],  S*num.yc[j] ]); params = params).p      for i in axes(num.Pt, 1), j in axes(num.Pt, 2)]

            push!(ε[:vx], mean(abs.(Vx_an .- num.Vx)))
            push!(ε[:vy], mean(abs.(Vy_an .- num.Vy)))
            push!(ε[:p],  mean(abs.(P_an  .- num.Pt)))
        end

        for (sym, label, color) in fields
            sc = scatter!(ax, invh, ε[sym]; color = color, markersize = 9)
            if k == 1
                push!(plots_for_legend, sc)
                push!(labels_for_legend, label)
            end
        end

        # O(1) reference line
        xs   = [minimum(invh), maximum(invh)]
        y0   = maximum(ε[sym][1] for (sym, _, _) in fields)
        C    = xs[1]*y0*1.3
        ln   = lines!(ax, xs, C ./ xs; color = :royalblue, linewidth = 2)
        if k == 1
            push!(plots_for_legend, ln)
            push!(labels_for_legend, "O(1)")
        end
    end

    Legend(fig[3, 1:2], plots_for_legend, labels_for_legend; orientation = :horizontal, nbanks = 1)

    save("convergence_JoukowskyDMvsNumerics_v2.png", fig, px_per_unit = 2)
    display(fig)
end
=#

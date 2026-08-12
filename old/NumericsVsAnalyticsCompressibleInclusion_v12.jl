using CairoMakie, MathTeXEngine, Printf, Statistics, LinearAlgebra, ExactFieldSolutions, StaticArrays

include("./Analytics_v12.jl")

@views av4_harm(A) = 1.0./( 0.25.*(1.0./A[1:end-1,1:end-1].+1.0./A[2:end,1:end-1].+1.0./A[1:end-1,2:end].+1.0./A[2:end,2:end]) ) 

# can be replaced by AD
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

# 2D Stokes routine
@views function Stokes2D(n, test)

    # Physics
    Lx, Ly   = 1.0, 1.0     # domain size
    ri       = 0.1          # inclusion radius
    comp     = true         # keep true :D     
    one_iter = false        # true: simple DR --- false: PH/DR

    if test == :Schmid2003
        # Schmid2003: pure shear, incompressible
        ηm, ηi = 1.0, 1e-2
        ξm, ξi = 1e10, 1e10
        ε̇ =-1e0
        γ̇ = 0.0
        ζ̇ = 0.0
    elseif test==:Duretz2026_1
        # Test 1: pure shear WITH finite compressibility
        ηm, ηi = 1.0, 1e-2
        ξm, ξi = 1e0, 1e0
        ε̇ =-1e0
        γ̇ = 0.0
        ζ̇ = 0.0
    elseif test==:Duretz2026_2
        # Test 2: pure shear + expansion WITH finite compressibility
        ηm, ηi = 1.0, 1e-1
        ξm, ξi = 1e0, 1e0
        ε̇ = 1e-10
        γ̇ = 1.0
        ζ̇ = 1e0
    elseif test==:Duretz2026_3
        # Test 3: pure shear + expansion + simple shear WITH finite compressibility
        ηm, ηi = 1.0, 1e-1
        ξm, ξi = 1e2, 1e2
        ε̇ = 1e0
        γ̇ = 1.0
        ζ̇ = 1e0
    end

    if comp
        f_anal = Analytics
        params = (ηm=1.0, ηi=ηi, ξm=ξm, ξi=ξm, ri=ri, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇)
    else
        # Call either the solution from ExactFieldSolutions
        # f_anal = Stokes2D_Schmid2003
        # params = (mm = 1.0, mc = ηi, rc =ri, gr = 0.0, er = ε̇)
        # ... or call the "new" solution with "infinite" bulk viscosities
        f_anal = Analytics
        params = (ηm=1.0, ηi=ηi, ξm=1e100*ξm, ξi=1e100*ξm, ri=ri, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇)
    end

    # Numerics
    ncx, ncy = n, n         # numerical grid resolution
    ϵ        = 1e-7         # tolerance
    iterMax  = 1e5          # max number of iters
    nout     = 100          # check frequency
    c_fact   = 0.5          # damping factor
    dτ_local = true         # helps a little bit sometimes, sometimes not! 
    γfact    = 60           # penalty: multiplier to the arithmetic mean of η
    rel_drop = 1e-3         # relative drop of velocity residual per PH iteration
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
    # Initialisation
    xce, yce = LinRange(-Lx/2-Δx/2, Lx/2+Δx/2, ncx+2), LinRange(-Ly/2-Δy/2, Ly/2+Δy/2, ncy+2)
    xc, yc   = LinRange(-Lx/2+Δx/2, Lx/2-Δx/2, ncx), LinRange(-Ly/2+Δy/2, Ly/2-Δy/2, ncy)
    xv, yv   = LinRange(-Lx/2, Lx/2, ncx+1), LinRange(-Ly/2, Ly/2, ncy+1)
    # Material heterogeneity: circular inclusion 
    ηv_sharp   .= ηm
    ηv_sharp[(xv).^2 .+ (yv').^2 .< ri^2 ] .= ηi
    ηc_sharp   .= ηm
    ηc_sharp[(xc).^2 .+ (yc').^2 .< ri^2 ] .= ηi
    # Use sharp viscosity
    ηc .= ηc_sharp
    ηv .= ηv_sharp
    # Harmonic averaging mimicking PIC interpolation
    ηc    .= av4_harm(ηv_sharp)
    ηv[2:end-1,2:end-1] .= av4_harm(ηc_sharp)

    # Bulk viscosity
    ηb    .= ξm

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
        dτVx =  2.0./sqrt.(λmaxVx)*0.99
        dτVy =  2.0./sqrt.(λmaxVy)*0.99
    else
        dτVx =  2.0./sqrt.(maximum(λmaxVx))*0.99 
        dτVy =  2.0./sqrt.(maximum(λmaxVy))*0.99
    end
    βVx .= 2 .* dτVx ./ (2 .+ cVx.*dτVx)
    βVy .= 2 .* dτVy ./ (2 .+ cVy.*dτVy)
    αVx .= (2 .- cVx.*dτVx) ./ (2 .+ cVx.*dτVx)
    αVy .= (2 .- cVy.*dτVy) ./ (2 .+ cVy.*dτVy)
    # Initial condition
    Vx     .=   ε̇.*xv .+  0*yce' .+  ζ̇ .*xv
    Vy     .=   0*xce .-  ε̇.*yv' .+  ζ̇ .*yv'
    Vx_ini  = copy(Vx)
    Vy_ini  = copy(Vy)
    Vx[2:end-1,:] .= 0 # ensure non zero initial pressure residual
    Vy[:,2:end-1] .= 0 # ensure non zero initial pressure residual
    Va = (
        x = zeros(ncx+1, ncy+2),
        y = zeros(ncx+2, ncy+1)
    )
    @time for I in CartesianIndices(Vx)
        i, j      = I[1], I[2]
        X         = @SVector([xv[i]; yce[j]])
        sol       = f_anal(X; params=params)
        Vx[i,j]   = sol.V[1] 
        Va.x[i,j] = sol.V[1] 
    end

    for I in CartesianIndices(Vy)
        i, j      = I[1], I[2]
        X         = @SVector([xce[i]; yv[j]])
        sol       = f_anal(X; params=params)
        Vy[i,j]   = sol.V[2] 
        Va.y[i,j] = sol.V[2] 
    end

    VxS, VxN = zeros(size(Vx,1)), zeros(size(Vx,1))
    for i in eachindex(VxS)
        X      = @SVector([xv[i]; -1/2])
        sol    = f_anal(X; params=params)
        VxS[i] = sol.V[1] 
        X      = @SVector([xv[i];  1/2])
        sol    = f_anal(X; params=params)
        VxN[i] = sol.V[1] 
    end

    VyW, VyE = zeros(size(Vy,2)), zeros(size(Vy,2))
    for j in eachindex(VyW)
        X      = @SVector([-1/2, yv[j]])
        sol    = f_anal(X; params=params)
        VyW[j] = sol.V[2] 
        X      = @SVector([1/2; yv[j]])
        sol    = f_anal(X; params=params)
        VyE[j] = sol.V[2] 
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
        X          = @SVector([xc[i]; yc[j]])
        sol        = f_anal(X; params=params)
        Pa[i,j]    = sol.p  
        Sa.xx[i,j] = sol.τ[1,1] - sol.p 
        Sa.yy[i,j] = sol.τ[2,2] - sol.p  
        Sa.zz[i,j] = sol.τzz    - sol.p  
    end

    for I in CartesianIndices(Sa.xy)
        i, j        = I[1], I[2]
        X           = @SVector([xv[i]; yv[j]])
        sol         = f_anal(X; params=params)
        Sa.xy[i,j]  = sol.τ[1,2] 
    end

    # Iteration loop
    errVx0 = 1.0;  errVy0 = 1.0;  errPt0 = 1.0 
    errVx00= 1.0;  errVy00= 1.0; 
    iter=1; err=2*ϵ; err_evo_V=[]; err_evo_P=[]; err_evo_it=[]
    @time for itPH = 1:50
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
    # Total stress
    Sxx = -Pt .+ Txx
    Syy = -Pt .+ Tyy
    Szz = -Pt .+ (-Txx .- Tyy)
    # Errors
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
    @show mean(Ve.x)
    @show mean(Ve.y)
    @show mean(Pe)
    @show mean(Se.xx)
    @show mean(Se.yy)
    @show mean(Se.zz)
    @show mean(Se.xy)
    # Plotting
    fig = Figure(size=(900, 1000))
    # -------------
    ax  = Axis(fig[1,1], aspect=DataAspect())
    hm  = heatmap!(ax, xv, yce, Vx , colormap=(Reverse(:matter), 1))
    ax  = Axis(fig[2,1], aspect=DataAspect())
    hm  = heatmap!(ax, xce, yv, Vy , colormap=(Reverse(:matter), 1))
    ax  = Axis(fig[3,1], aspect=DataAspect())
    hm  = heatmap!(ax, xc, yc, Pt, colormap=(Reverse(:matter), 1))
    ax  = Axis(fig[4,1], aspect=DataAspect())
    hm  = heatmap!(ax, xc, yc, Sxx, colormap=(Reverse(:matter), 1))
    ax  = Axis(fig[5,1], aspect=DataAspect())
    hm  = heatmap!(ax, xc, yc, Szz, colormap=(Reverse(:matter), 1))
    ax  = Axis(fig[6,1], aspect=DataAspect())
    hm  = heatmap!(ax, xv, yv, Txy, colormap=(Reverse(:matter), 1))
    # -------------
    ax  = Axis(fig[1,2], aspect=DataAspect())
    hm  = heatmap!(ax, xv, yce, Va.x , colormap=(Reverse(:matter), 1))
    Colorbar(fig[1,3], hm)
    ax  = Axis(fig[2,2], aspect=DataAspect())
    hm  = heatmap!(ax, xce, yv, Va.y , colormap=(Reverse(:matter), 1))
    Colorbar(fig[2,3], hm)
    ax  = Axis(fig[3,2], aspect=DataAspect())
    hm  = heatmap!(ax, xc, yc, Pa, colormap=(Reverse(:matter), 1))
    Colorbar(fig[3,3], hm)
    ax  = Axis(fig[4,2], aspect=DataAspect())
    hm  = heatmap!(ax, xc, yc, Sa.xx, colormap=(Reverse(:matter), 1))
    Colorbar(fig[4,3], hm)
    ax  = Axis(fig[5,2], aspect=DataAspect())
    hm  = heatmap!(ax, xc, yc, Sa.zz, colormap=(Reverse(:matter), 1))
    Colorbar(fig[5,3], hm)
    ax  = Axis(fig[6,2], aspect=DataAspect())
    hm  = heatmap!(ax, xv, yv, Sa.xy, colormap=(Reverse(:matter), 1))
    Colorbar(fig[6,3], hm)
    # -------------
    ax  = Axis(fig[1,4], aspect=DataAspect())
    hm  = heatmap!(ax, xv, yce, Ve.x, colormap=(Reverse(:matter), 1))
    Colorbar(fig[1,5], hm)
    ax  = Axis(fig[2,4], aspect=DataAspect())
    hm  = heatmap!(ax, xce, yv, Ve.y, colormap=(Reverse(:matter), 1))
    Colorbar(fig[2,5], hm)
    ax  = Axis(fig[3,4], aspect=DataAspect())
    hm  = heatmap!(ax, xc, yc, Pe, colormap=(Reverse(:matter), 1))
    Colorbar(fig[3,5], hm)
    ax  = Axis(fig[4,4], aspect=DataAspect())
    hm  = heatmap!(ax, xc, yc, Se.xx, colormap=(Reverse(:matter), 1))
    Colorbar(fig[4,5], hm)
    ax  = Axis(fig[5,4], aspect=DataAspect())
    hm  = heatmap!(ax, xc, yc, Se.zz, colormap=(Reverse(:matter), 1))
    Colorbar(fig[5,5], hm)
    ax  = Axis(fig[6,4], aspect=DataAspect())
    hm  = heatmap!(ax, xv, yv, Se.xy, colormap=(Reverse(:matter), 1))
    Colorbar(fig[6,5], hm)
    # -------------
    display(fig)
    @show iter/ncx
    L1 = (Vx=mean(Ve.x), Vy=mean(Ve.y), P=mean(Pe), σxx=mean(Se.xx), σyy=mean(Se.yy), σzz=mean(Se.zz), σxy=mean(Se.xy))
    return L1
end

# Stokes2D( 51, :Duretz2026_3 )

let 

     #test = :Schmid2003
     #test = :Duretz2026_1   
     test = :Duretz2026_2
    #test = :Duretz2026_3

    nc = [200 ]
    L1  = (
        Vx  = zeros(length(nc)), 
        Vy  = zeros(length(nc)), 
        P   = zeros(length(nc)), 
        σxx = zeros(length(nc)), 
        σyy = zeros(length(nc)), 
        σzz = zeros(length(nc)), 
        σxy = zeros(length(nc))
    )
    for i in eachindex(nc) 
        errors = Stokes2D( nc[i], test )
        L1.Vx[i]  = errors.Vx  
        L1.Vy[i]  = errors.Vy  
        L1.P[i]   = errors.P   
        L1.σxx[i] = errors.σxx 
        L1.σyy[i] = errors.σyy 
        L1.σzz[i] = errors.σzz 
        L1.σxy[i] = errors.σxy
    end

    @show L1.Vx[1:end-1]  ./ L1.Vx[2:end] 
    @show L1.Vy[1:end-1]  ./ L1.Vy[2:end] 
    @show L1.P[1:end-1]   ./ L1.P[2:end] 
    @show L1.σxx[1:end-1] ./ L1.σxx[2:end] 
    @show L1.σyy[1:end-1] ./ L1.σyy[2:end] 
    @show L1.σzz[1:end-1] ./ L1.σzz[2:end] 
    @show L1.σxy[1:end-1] ./ L1.σxy[2:end]
end

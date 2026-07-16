# Numerical pseudotransient solver
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

# Stokes solver
@views function Stokes2D(n, test; formulation=:new)

    f_anal = formulation == :new ? Analytics_new : Analytics_old

    Lx, Ly   = 1.0, 1.0
    comp     = true
    one_iter = false

    params = preset_circle(test)
    ηm, ηi, ξm, ξi, γ̇, ε̇, ζ̇, ri = params

    # For BCs use the selected analytical solution with physical compressibility
<<<<<<< Updated upstream:src/Numerics.jl
    bc_params = comp ? (ηm=ηm, ηi=ηi, ξm=ξm, ξi=ξm, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇, ri=ri) :
                       (ηm=ηm, ηi=ηi, ξm=1e100*ξm, ξi=1e100*ξm, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇, ri=ri)
=======
    bc_params = comp ? (ηm=ηm, ηi=ηi, ξm=ξm, ξi=ξi, ri=ri, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇) :
                   (ηm=ηm, ηi=ηi, ξm=1e100*ξm, ξi=1e100*ξi, ri=ri, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇)
>>>>>>> Stashed changes:NumericsVsAnalytics_systematics.jl

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
    ηb[(xc).^2 .+ (yc').^2 .< ri^2] .= ξi
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
    Sa = (
        xx = zeros(ncx, ncy),
        yy = zeros(ncx, ncy),
        zz = zeros(ncx, ncy),
        xy = zeros(ncx+1, ncy+1),
    )
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
    X  = (xv=xv, yv=yv, xce=xce, yce=yce, xc=xc, yc=yc)

    @show iter/ncx
    return L1, V, Pt, S, Va, Pa, Sa, X
end
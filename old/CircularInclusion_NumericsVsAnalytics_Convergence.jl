using InclusionAnalytics
using CairoMakie, MathTeXEngine
using ExactFieldSolutions, JLD2

# Compare solutions
let
    formulation = :new             #:old | :new
    Lx, Ly = 1.0, 1.0
    nc          = [51, 101, 201]#, 401]              #resolution

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

        errors, V, P, S, Va, Pa, Sa, X =
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

    for itest in eachindex(tests)
        @info "Error reduction: $(itest)"
        @show L1.Vx[itest, 1:end-1] ./ L1.Vx[itest, 2:end] 
    end
    # jldsave("TruncationSystematics.jld2", errors=L1)
end

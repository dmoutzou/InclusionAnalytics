function eval_analytics_stag(f_anal, X, params, geometry, scale)
    # Analytical reference on the same staggered grids as the numerics
    ncx, ncy = length(X.xc), length(X.xc)
    Va = (
        x = zeros(ncx+1, ncy+2),
        y = zeros(ncx+2, ncy+1),
    )
    Pa = zeros(ncx+0, ncy+0)
    τa = (
        xx = zeros(ncx+0, ncy+0),
        yy = zeros(ncx+0, ncy+0),
        zz = zeros(ncx+0, ncy+0),
        xy = zeros(ncx+1, ncy+1),
    )
    for I in CartesianIndices(Va.x)
        i, j = I[1], I[2]
        Va.x[i,j] = f_anal(@SVector([X.xv[i]; X.yce[j]]); params=params, geometry=geometry, S=scale).V[1]
    end
    for I in CartesianIndices(Va.y)
        i, j = I[1], I[2]
        Va.y[i,j] = f_anal(@SVector([X.xce[i]; X.yv[j]]); params=params, geometry=geometry, S=scale).V[2]
    end
    for I in CartesianIndices(Pa)
        i, j = I[1], I[2]
        sol = f_anal(@SVector([X.xc[i]; X.yc[j]]); params=params, geometry=geometry, S=scale)
        Pa[i,j]    = sol.p
        τa.xx[i,j] = sol.τ[1,1] - sol.p
        τa.yy[i,j] = sol.τ[2,2] - sol.p
        τa.zz[i,j] = sol.τzz    - sol.p
    end
    for I in CartesianIndices(τa.xy)
        i, j = I[1], I[2]
        τa.xy[i,j] = f_anal(@SVector([X.xv[i]; X.yv[j]]); params=params, geometry=geometry, S=scale).τ[1,2]
    end

    return Va, Pa, τa 
end
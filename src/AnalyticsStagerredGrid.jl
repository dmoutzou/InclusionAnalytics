function analytics_stag(f_anal, X, params, geometry)
    # Analytical reference on the same staggered grids as the numerics
    ncx, ncy = length(X.xc), length(X.xc)
    Va = (
        x = zeros(ncx+1, ncy+2),
        y = zeros(ncx+2, ncy+1),
    )
    Pa = zeros(ncx+0, ncy+0)
    σa = (
        xx = zeros(ncx+0, ncy+0),
        yy = zeros(ncx+0, ncy+0),
        zz = zeros(ncx+0, ncy+0),
        xy = zeros(ncx+1, ncy+1),
    )
    for I in CartesianIndices(Va.x)
        i, j = I[1], I[2]
        Va.x[i,j] = f_anal(@SVector([X.xv[i]; X.yce[j]]), params; geometry=geometry).V[1]
    end
    for I in CartesianIndices(Va.y)
        i, j = I[1], I[2]
        Va.y[i,j] = f_anal(@SVector([X.xce[i]; X.yv[j]]), params; geometry=geometry).V[2]
    end
    for I in CartesianIndices(Pa)
        i, j = I[1], I[2]
        sol = f_anal(@SVector([X.xc[i]; X.yc[j]]), params; geometry=geometry)
        Pa[i,j]    = sol.p
        σa.xx[i,j] = sol.τ[1,1] - sol.p
        σa.yy[i,j] = sol.τ[2,2] - sol.p
        σa.zz[i,j] = sol.τzz    - sol.p
    end
    for I in CartesianIndices(σa.xy)
        i, j = I[1], I[2]
        σa.xy[i,j] = f_anal(@SVector([X.xv[i]; X.yv[j]]), params; geometry=geometry).τ[1,2]
    end

    return Va, Pa, σa 
end

function errors_stag(V, P, σ, Va, Pa, σa)
    # Errors
    Ve  = (x=abs.(V.x .- Va.x), y = abs.(V.y .- Va.y))
    Pe  = abs.(P .- Pa)
    σe = ( xx=abs.(σ.xx .- σa.xx), yy = abs.(σ.yy .- σa.yy), zz = abs.(σ.zz .- σa.zz),  xy = abs.(σ.xy .- σa.xy))
    @printf("mean|Vx error| = %1.3e   mean|Vy error| = %1.3e   mean|P error| = %1.3e\n", mean(Ve.x), mean(Ve.y), mean(Pe))
    @printf("mean|σxx error| = %1.3e   mean|σyy error| = %1.3e   mean|σzz error| = %1.3e   mean|σxy error| = %1.3e\n",
            mean(σe.xx), mean(σe.yy), mean(σe.zz), mean(σe.xy))
    return Ve, Pe, σe
end

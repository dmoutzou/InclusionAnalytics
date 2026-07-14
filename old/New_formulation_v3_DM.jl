using CairoMakie

#Physical parameters
ηm = 1.0;  ηi = 0.01
ξm = 1e-2;  ξi = 1e-2
ri = 0.1;  γ̇ = 0.0;  ε̇ = -1.0;  ζ̇ = 0.0
#Kolosov constants
νm = (3ξm - 2ηm) / (2*(3ξm + ηm))
νi = (3ξi - 2ηi) / (2*(3ξi + ηi))
κm = 3 - 4νm
κi = 3 - 4νi
@show νm, νi, κm, κi
#Far-field potentials (coming from boundary conditions)
P_m = -im*ηm*γ̇/(κm+1) + 2*ηm*ζ̇/(κm-1) #BCs related to φ potentials
Q_m = (im*γ̇ - 2ε̇)*ηm                  #BCs related to ψ potentials
#Take the real and imaginary part of P_m (will be needed)
P_mR = real(P_m)
P_mI = imag(P_m)
@show P_m, Q_m
#S1: B1R from real part of V1
B1R = (κm+1)*ηi*P_mR / (2ηi + (κi-1)*ηm)
#S2: B1I from imaginary part of V1
B1I = (κm+1)*ηi*P_mI / ((κi+1)*ηm)
B1  = B1R + im*B1I
#S3: A1 from e^{-iθ} velocity condition
A1  = (ηi - ηm)*conj(Q_m)*ri^2 / (κm*ηi + ηm)
#T3: A3 from e^{3iθ} traction condition
A3  = ri^2 * A1
#T2: B2 from e^{-iθ} traction condition
B2  = conj(A1/ri^2 + conj(Q_m))
#S4: A2 from T1 substituting S1
A2  = 2ri^2 * ((κm-1)*ηi - (κi-1)*ηm) / (2ηi + (κi-1)*ηm) * P_mR
@show B1, A1, A2, A3, B2
#Field evaluation
function Analytical_Solution_new(x, y)
    z   = x + im*y
    z̄   = conj(z)
    r   = abs(z)
    if r > ri
        ϕ   =  P_m*z  + A1/z
        ϕ′  =  P_m    - A1/z^2
        ϕ′′ =          2A1/z^3
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
    #Get stresses
    Qval = z̄*ϕ′′ + ψ′
    sxx  = 2real(ϕ′) - real(Qval)
    syy  = 2real(ϕ′) + real(Qval)
    sxy  = imag(Qval)
    szz  = ν_loc*(sxx + syy)
    p    = -(sxx + syy + szz) / 3
    #Velocity calculation
    vel  = (κ_loc*ϕ - z*conj(ϕ′) - conj(ψ)) / (2*η_loc)
    return p, real(vel), imag(vel), sxx, syy, sxy
end
#Spatial grid
n  = 300
x  = LinRange(-0.5, 0.5, n)
y  = LinRange(-0.5, 0.5, n)
P  = zeros(n,n)
Vx = zeros(n,n)
Vy = zeros(n,n)
for (j,yj) in enumerate(y), (i,xi) in enumerate(x)
    P[i,j], Vx[i,j], Vy[i,j] = Analytical_Solution_new(xi, yj)[1:3]
end
#print some points
println("\nField values at selected points:")
println("Maximum Vx is ", maximum(Vx))
println("Inside  (0,    0):    p, vx, vy = ", Analytical_Solution_new(0.0,  0.0)[1:3])
println("On rim  (ri,   0):    p, vx, vy = ", Analytical_Solution_new(ri*1.001, 0.0)[1:3])
println("Outside (2ri,  0):    p, vx, vy = ", Analytical_Solution_new(2ri, 0.0)[1:3])
println("Outside (0,   2ri):   p, vx, vy = ", Analytical_Solution_new(0.0, 2ri)[1:3])
#Plotting
clim = maximum(abs.(P))
vlim = maximum(abs.(Vx))
fig  = Figure(size=(1400, 500))
ax1  = Axis(fig[1,1], aspect=DataAspect(), title="Pressure")
hm1  = heatmap!(ax1, x, y, P, colormap=(Reverse(:matter), 1))
Colorbar(fig[1,2], hm1, label="p")
ax2  = Axis(fig[1,3], aspect=DataAspect(), title="Vx")
hm2  = heatmap!(ax2, x, y, Vx, colormap=(Reverse(:matter), 1))
Colorbar(fig[1,4], hm2, label="vx")
ax3  = Axis(fig[1,5], aspect=DataAspect(), title="Vy")
hm3  = heatmap!(ax3, x, y, Vy, colormap=(Reverse(:matter), 1))
Colorbar(fig[1,6], hm3, label="vy")
display(fig)
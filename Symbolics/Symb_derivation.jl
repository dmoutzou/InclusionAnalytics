# Symbolic derivation of the mode-matching equations and closed-form coefficients
# for a compressible host-inclusion problem. 
#
# Ansatz 
#
#   Host:
#       phi_m(zeta) = Pm*omega(zeta) + A1/zeta
#       psi_m(zeta) = Qm*omega(zeta) + A2/zeta + A3/(zeta^3 - zeta)
#
#   Inclusion:
#       phi_i(zeta) = B1*omega(zeta)
#       psi_i(zeta) = B2*omega(zeta)
#
#   Traction:  phi + [omega/conj(omega')]*conj(phi') + conj(psi)   continuous
#   Velocity:  ( kappa*phi - [omega/conj(omega')]*conj(phi') - conj(psi) ) / eta   continuous
#
# NOTE: the chosen ansatz produces a 12x10 system (12 equations but 10 unknowns). However, 
# 2 of these equations are dependent and can be linearly be produced by the other 10. Also the "solve"
# function in SYmbolics.jl only accepts square matrices. So I reduce the system by finding which equations
# do not contribute to the solution numerically

using Symbolics, LinearAlgebra, Random
Random.seed!(1)   

# 1. Variables

@variables zeta::Real r::Real
@variables etam::Real etai::Real kappam::Real kappai::Real
@variables Pm Qm Pmc Qmc A1 A2 A3 B1 B2 A1c A2c A3c B1c B2c

#Real/imaginary decompositions 
@variables pmx::Real pmy::Real qmx::Real qmy::Real
@variables a1x::Real a1y::Real a2x::Real a2y::Real a3x::Real a3y::Real b1x::Real b1y::Real b2x::Real b2y::Real

# 2. Potentials and derivatives
omega = zeta + 1/zeta

phi_m = Pm*omega + A1/zeta
psi_m = Qm*omega + A2/zeta + A3/(zeta^3 - zeta)
phi_i = B1*omega
psi_i = B2*omega

phi_m_p = Symbolics.derivative(phi_m, zeta)
phi_i_p = Symbolics.derivative(phi_i, zeta)

#3. Boundary conjugates
conjSubsM = Dict(Pm => Pmc, Qm => Qmc, A1 => A1c, A2 => A2c, A3 => A3c)
conjSubsI = Dict(B1 => B1c, B2 => B2c)

boundary_conj_m(expr) = substitute(substitute(expr, conjSubsM), Dict(zeta => r^2/zeta))
boundary_conj_i(expr) = substitute(substitute(expr, conjSubsI), Dict(zeta => r^2/zeta))

phi_m_p_conj = boundary_conj_m(phi_m_p)
psi_m_conj   = boundary_conj_m(psi_m)
phi_i_p_conj = boundary_conj_i(phi_i_p)
psi_i_conj   = boundary_conj_i(psi_i)

omega_pbar = 1 - zeta^2/r^4      # = conj(omega'(zeta)) on |zeta| = r
Omega      = omega/omega_pbar

# 4. Traction and velocity continuity
traction_residual = (phi_m + Omega*phi_m_p_conj + psi_m_conj) -
                     (phi_i + Omega*phi_i_p_conj + psi_i_conj)

velocity_residual = (kappam*phi_m - Omega*phi_m_p_conj - psi_m_conj)/etam -
                     (kappai*phi_i - Omega*phi_i_p_conj - psi_i_conj)/etai

#Clear the conj(omega') 
traction_cleared = expand(simplify(traction_residual * omega_pbar * zeta^3))
velocity_cleared = expand(simplify(velocity_residual * omega_pbar * etam * etai * zeta^3))

# 5. Extract all the modes
trac_dict, trac_rem = Symbolics.polynomial_coeffs(traction_cleared, (zeta,))
vel_dict,  vel_rem  = Symbolics.polynomial_coeffs(velocity_cleared, (zeta,))

println("Traction: harmonics present (powers of zeta): ", collect(keys(trac_dict)))
println("Traction remainder (should be 0): ", trac_rem)
println("Velocity: harmonics present (powers of zeta): ", collect(keys(vel_dict)))
println("Velocity remainder (should be 0): ", vel_rem)

# 6. Build the full system of equations
complexSyms = [A1, A1c, A2, A2c, A3, A3c, B1, B1c, B2, B2c, Pm, Pmc, Qm, Qmc]
complexVals = [a1x+im*a1y, a1x-im*a1y, a2x+im*a2y, a2x-im*a2y, a3x+im*a3y, a3x-im*a3y,
               b1x+im*b1y, b1x-im*b1y, b2x+im*b2y, b2x-im*b2y,
               pmx+im*pmy, pmx-im*pmy, qmx+im*qmy, qmx-im*qmy]

unknowns = [a1x, a1y, a2x, a2y, a3x, a3y, b1x, b1y, b2x, b2y]
params   = [r, etam, etai, kappam, kappai]

allCoeffs = vcat(collect(values(trac_dict)), collect(values(vel_dict)))
fc  = eval(build_function(allCoeffs, complexSyms, params)[1])
cks = Base.invokelatest(fc, complexVals, params)

eqs = Num[]
for ck in cks
    rePart = simplify(expand(real(ck)))
    imPart = simplify(expand(imag(ck)))
    iszero(rePart) || push!(eqs, rePart)
    iszero(imPart) || push!(eqs, imPart)
end

println("\nTotal nontrivial real equations extracted: ", length(eqs))
println("Number of real unknowns: ", length(unknowns))

# 7. Consistency check to find the dependent equations
M = Symbolics.jacobian(eqs, unknowns)
subst  = Dict(p => 0.5 + rand()*3.5 for p in params)   # generic point, away from 0/degeneracies
Mnum   = Float64.(Symbolics.value.(substitute.(M, (subst,))))
rk     = rank(Mnum)

println("\nCoefficient matrix size: ", size(M), "   numeric rank: ", rk)
if rk == length(unknowns)
    println("=> System is CONSISTENT and uniquely determined for general kappam, kappai.")
else
    error("=> WARNING: system is under/over-determined -- ansatz needs revisiting.")
end

pivotRows = sort(qr(Matrix(Mnum'), ColumnNorm()).p[1:rk])
eqs_sq    = eqs[pivotRows]

#   8. Solve in closed form
sol = Symbolics.symbolic_linear_solve(eqs_sq .~ 0, unknowns)
d   = Dict(unknowns .=> sol)

A1_sol = d[a1x] + im*d[a1y]
A2_sol = d[a2x] + im*d[a2y]
A3_sol = d[a3x] + im*d[a3y]
B1_sol = d[b1x] + im*d[b1y]
B2_sol = d[b2x] + im*d[b2y]

# 9. Here I make sure that just the simplifications I was doing are correct
r4 = r^4

K   = -2etai*kappam*pmx + 2etai*pmx + 2etai*qmx*r^2 + 2etam*kappai*pmx +
       etam*kappai*qmx*r^2 - 2etam*pmx - etam*qmx*r^2
L   = etai*kappam*pmx*r4 - etai*pmx - etai*qmx*r^2 + etam*pmx*r4 + etam*pmx + etam*qmx*r^2
Mm  = etai*kappam - etai - etam*kappai + etam
DEN = 2etai^2*kappam*r4 - 2etai^2*kappam + etai*etam*kappai*kappam*r4 + etai*etam*kappai -
      etai*etam*kappam*r4 + 2etai*etam*kappam + 2etai*etam*r4 - etai*etam +
      etam^2*kappai*r4 - etam^2*kappai - etam^2*r4 + etam^2
den1 = etai*kappam*r4 + etai + etam*r4 - etam
den2 = etai*kappam*r4 - etai + etam*r4 + etam

A1 = (etai - etam)*(r4 - 1)*K/DEN + im*( qmy*r^2*(etam - etai)*(r4 - 1)/den1 )
A2 = 2*(r4 - 1)*Mm*L/(r^2*DEN) + im*0
A3 = (etai - etam)*(r4 - 1)*(r4 + 1)*K/(r^2*DEN) + im*( qmy*(etam - etai)*(r^8 - 1)/den1 )
B1 = etai*(kappam + 1)*L/DEN +
     im*( etai*(kappam + 1)*(etai*kappam*pmy*r4 + etai*pmy + etai*qmy*r^2 +
          etam*pmy*r4 - etam*pmy - etam*qmy*r^2) / (etam*(kappai + 1)*den1) )
B2 = etai*r^2*(kappam + 1)*K/DEN + im*( etai*qmy*r4*(kappam + 1)/den1 )

function isnumzero(expr; ntrials=5, tol=1e-6)
    free = collect(Symbolics.get_variables(expr))
    isempty(free) && return abs(ComplexF64(Symbolics.value(expr))) < tol
    f = eval(build_function(expr, free))
    for _ in 1:ntrials
        v = Base.invokelatest(f, [0.2 + rand()*3 for _ in free])
        abs(ComplexF64(v)) > tol && return false
    end
    return true
end

println("\n--- Closed-form coefficients (factored form, matches analytics_ellipse) ---")
for (name, raw, compact) in (("A1", A1_sol, A1), ("A2", A2_sol, A2), ("A3", A3_sol, A3),
                              ("B1", B1_sol, B1), ("B2", B2_sol, B2))
    ok = isnumzero(expand(raw - compact))
    println(name, " = ", compact, "   [", ok ? "verified against raw solve" : "MISMATCH -- do not trust", "]")
end

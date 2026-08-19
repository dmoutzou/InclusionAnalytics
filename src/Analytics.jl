# physical position from the conformal coordinate
joukowski(λ) = λ + 1.0/λ          # ω(λ) = z   (change here if your map differs)
t_to_ri(t) = sqrt((t - 1.0)*(t + 1.0)) / (t - 1.0)

# The numerical solver is using the true cartesian coordinates, so I need the inverse mapping of the Joukowsky transform
function inv_joukowski(z)
    disc = sqrt(z^2 - 4.0 + 0im)
    ζ1   = (z + disc)/2
    ζ2   = (z - disc)/2
    return abs(ζ1) >= abs(ζ2) ? ζ1 : ζ2
end
to_zeta(X) = (ζ = inv_joukowski(complex(X[1], X[2])); @SVector([real(ζ), imag(ζ)]))

# Function to acquire the axes of the inclusion based on t (
function ellipse_axes(t)
    ri = t_to_ri(t)
    a  = ri + 1.0/ri
    b  = ri - 1.0/ri
    return a, b
end


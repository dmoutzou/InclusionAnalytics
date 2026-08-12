function preset(test, geometry)
    # if geometry == :elliptical
        if test == :Schmid2003
            ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e10, 1e10
            ε̇ = -0.5;  γ̇ = 0.0;  ζ̇ = 0.0;  ε̇zz = 0.0
            ri = 0.1;  t = 2.0;   α = 0.0;
        elseif test == :Duretz2026_1_ps
            ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e0, 10e0
            ε̇ = -1e0;  γ̇ = 0.0;  ζ̇ = 0.0;  ε̇zz = 0.0
            ri = 0.1;  t = 2.0;   α = 0.0;
        elseif test == :Duretz2026_2_ss
            ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e0, 10e0
            ε̇ = 0.0; γ̇ = -1.0;  ζ̇ = 0.0;  ε̇zz = 0.0
            ri = 0.1;  t = 2.0;   α = 0.0;
        elseif test == :Duretz2026_3_exp
            ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e0, 10e0
            ε̇ = 0.0;   γ̇ = 0.0;  ζ̇ = 1.0;  ε̇zz = 0.0
            ri = 0.1;  t = 2.0;   α = 0.0;
        elseif test == :Duretz2026_4_comp
            ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e0, 10e0
            ε̇ = 0.0;   γ̇ = 0.0;  ζ̇ = -1.0;  ε̇zz = 0.0
            ri = 0.1;  t = 2.0;   α = 0.0;
        elseif test == :Duretz2026_5_mixed
            ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e0, 10e0
            ε̇ = 1.0;   γ̇ = 1.0;  ζ̇ = 2.0;  ε̇zz = 0.0
            ri = 0.1;  t = 2.0;   α = 0.0;
        elseif test == :Duretz2026_6_angled
            ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e0, 10e0
            ε̇ = 0.0;   γ̇ = 1.0;  ζ̇ = 0.0;  ε̇zz = 0.0
            ri = 0.1;  t = 2.0;   α = 30/180*pi;
        elseif test == :Duretz2026_7_oop
            ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e0, 10e0
            ε̇ = 0.0;   γ̇ = 0.0;  ζ̇ = 1.0;  ε̇zz = 1.0
            ri = 0.1;  t = 2.0;   α = 0.0;
        elseif test == :Duretz2026_8_zeroop
            ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e0, 10e0
            ε̇ = 0.0;   γ̇ = 0.0;  ζ̇ = 1.0;  ε̇zz = 0.0
            ri = 0.1;  t = 2.0;   α = 0.0;
        else
            error("unknown test $test")
        end
        # return (ηm=ηm, ηi=ηi, ξm=ξm, ξi=ξi, t=t, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇, α=α)

    # elseif geometry == :circular
    #     if test == :Schmid2003
    #         ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e10, 1e10
    #         ε̇ = -1e0;  γ̇ = 0.0;  ζ̇ = 0.0
    #         ri = 0.1
    #     elseif test == :Duretz2026_1
    #         ηm, ηi = 1.0, 1e-2;  ξm, ξi = 1e0, 1e0
    #         ε̇ = -1e0;  γ̇ = 0.0;  ζ̇ = 0.0
    #         ri = 0.1
    #     elseif test == :Duretz2026_2
    #         ηm, ηi = 1.0, 1e-1;  ξm, ξi = 1e0, 1e0
    #         ε̇ = 1e-10; γ̇ = 1.0;  ζ̇ = 1e0
    #         ri = 0.1
    #     elseif test == :Duretz2026_3
    #         ηm, ηi = 1.0, 1e-1;  ξm, ξi = 1.0, 1.0
    #         ε̇ = 1e-10; γ̇ = 0.0;  ζ̇ = 1e0
    #         ri = 0.1
    #     else
    #         error("Unknown test: $test")
    #     end
        # return (ηm=ηm, ηi=ηi, ξm=ξm, ξi=ξi, γ̇=γ̇, ε̇=ε̇, ζ̇=ζ̇, ri=ri)

    # else
    #     error("Unknown geometry: $geometry")
    # end

    if geometry == :elliptical
        r1, r2     = ellipse_axes(t) # true physical semi-axes (a >= 2 always)
        sc         = r2 / ri 
    else
        r1, r2, sc = ri, ri, 1.0
    end

    return (ηm=ηm, ηi=ηi, ξm=ξm, ξi=ξi, r1=r1/sc, r2=r2/sc, t=t, α=α, sc=sc, ε̇=ε̇, γ̇=γ̇, ζ̇=ζ̇, ε̇zz=ε̇zz)
end
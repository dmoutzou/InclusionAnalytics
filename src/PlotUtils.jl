# number -> LaTeX fragment: 0.01 -> "10^{-2}", 250 -> "2.5\times10^{2}"
function fmtnum(x)
    x == 0 && return "0"
    e = floor(Int, log10(abs(x)))
    if -2 ≤ e ≤ 3
        v = round(x, sigdigits=3)
        return string(v == round(v) ? Int(round(v)) : v)
    else
        m = round(x / 10.0^e, sigdigits=3)
        return m == 1.0 ? "10^{$e}" : "$(m)\\times 10^{$e}"
    end
end

# Gtk4 GUI for InclusionAnalytics
#
#   include("interface/gtk_interface.jl")
#
# The window freezes (does not repaint / respond) while "Run solver" is
# computing the numerical solution.

import Pkg
Pkg.activate(@__DIR__)

using Gtk4, Gtk4.GLib
using CairoMakie
using Statistics, Printf

# ---------------------------------------------------------------------------
# Load the existing InclusionAnalytics package code directly 
# ---------------------------------------------------------------------------
const PKG_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(PKG_ROOT, "src", "InclusionAnalytics.jl"))
using .InclusionAnalytics

CairoMakie.activate!()
const MAKIE_SCREEN_CONFIG = CairoMakie.ScreenConfig(1.0, 1.0, :good, true, false, nothing)

function build_params(geometry::Symbol, ηm, ηi, ξm, ξi, ε̇, γ̇, ζ̇, ε̇zz, ri, t, α)
    if geometry == :elliptical
        r1, r2 = InclusionAnalytics.ellipse_axes(t)
        sc = r2 / ri
    else
        r1, r2, sc = ri, ri, 1.0
    end
    return (ηm=ηm, ηi=ηi, ξm=ξm, ξi=ξi, r1=r1/sc, r2=r2/sc, t=t, α=α, sc=sc,
            ε̇=ε̇, γ̇=γ̇, ζ̇=ζ̇, ε̇zz=ε̇zz)
end

f_anal_for(geometry::Symbol) = geometry == :circular ? InclusionAnalytics.analytics_circle :
                                                          InclusionAnalytics.analytics_ellipse

active_text(combo) = Gtk4.G_.get_active_text(combo)

function flush_gtk_events!()
    while ccall((:g_main_context_iteration, Gtk4.GLib.libglib), Cint, (Ptr{Cvoid}, Cint), C_NULL, false) != 0
    end
end

# ---------------------------------------------------------------------------
# Solve results state (kept so the field selector can redraw without re-solving)
# ---------------------------------------------------------------------------
mutable struct SolveState
    ready::Bool
    geometry::Symbol
    nc::Int
    X::Any
    V::Any; P::Any; σ::Any
    Va::Any; Pa::Any; σa::Any
    Ve::Any; Pe::Any; σe::Any
end
const state = SolveState(false, :elliptical, 0, nothing, nothing, nothing, nothing,
                          nothing, nothing, nothing, nothing, nothing, nothing)

function field_arrays(key::String)
    gx, gy, an, num, label, share = _raw_field_arrays(key)
    scale = 2 / (state.X.xv[end] - state.X.xv[1])
    return scale .* gx, scale .* gy, an, num, label, share
end

function _raw_field_arrays(key::String)
    X, P, Pa, σ, σa = state.X, state.P, state.Pa, state.σ, state.σa
    if key == "Vx"
        return X.xv,  X.yce, state.Va.x, state.V.x, "Vx", false
    elseif key == "Vy"
        return X.xce, X.yv,  state.Va.y, state.V.y, "Vy", false
    elseif key == "P"
        return X.xc, X.yc, Pa, P, "P", true
    elseif key == "sigma_xx"
        return X.xc, X.yc, σa.xx, σ.xx, "sigma_xx (total)", true
    elseif key == "sigma_yy"
        return X.xc, X.yc, σa.yy, σ.yy, "sigma_yy (total)", true
    elseif key == "sigma_zz"
        return X.xc, X.yc, σa.zz, σ.zz, "sigma_zz (total)", true
    elseif key == "sigma_xy"
        return X.xv, X.yv, σa.xy, σ.xy, "sigma_xy (total)", true
    elseif key == "tau_xx"
        return X.xc, X.yc, σa.xx .+ Pa, σ.xx .+ P, "tau_xx (deviatoric)", true
    elseif key == "tau_yy"
        return X.xc, X.yc, σa.yy .+ Pa, σ.yy .+ P, "tau_yy (deviatoric)", true
    elseif key == "tau_zz"
        return X.xc, X.yc, σa.zz .+ Pa, σ.zz .+ P, "tau_zz (deviatoric)", true
    elseif key == "tau_xy"
        return X.xv, X.yv, σa.xy, σ.xy, "tau_xy (deviatoric)", true
    elseif key == "tau_II"
        τa_xx, τ_xx = σa.xx .+ Pa, σ.xx .+ P
        τa_yy, τ_yy = σa.yy .+ Pa, σ.yy .+ P
        τa_zz, τ_zz = σa.zz .+ Pa, σ.zz .+ P
        τa_xy, τ_xy = σa.xy, σ.xy
        avg(A) = 0.25 .* (A[1:end-1,1:end-1] .+ A[1:end-1,2:end] .+ A[2:end,1:end-1] .+ A[2:end,2:end])
        τaII = sqrt.(0.5 .* (τa_xx.^2 .+ τa_yy.^2 .+ τa_zz.^2) .+ avg(τa_xy).^2)
        τII  = sqrt.(0.5 .* (τ_xx.^2  .+ τ_yy.^2  .+ τ_zz.^2 ) .+ avg(τ_xy).^2)
        return X.xc, X.yc, τaII, τII, "tau_II", true
    else
        error("unknown field $key")
    end
end

function build_figure(key::String)
    gx, gy, an, num, label, share = field_arrays(key)
    diff = abs.(num .- an)
    mag  = max(maximum(abs, an), maximum(abs, num), 1.0)
    crange = if share
        lo, hi  = extrema(vcat(vec(an), vec(num)))
        minspan = 0.02 * mag
        if hi - lo < minspan
            c = 0.5 * (lo + hi)
            lo, hi = c - 0.5 * minspan, c + 0.5 * minspan
        end
        (colorrange = (lo, hi),)
    else
        NamedTuple()
    end
    dmax = max(maximum(abs, diff), 0.01 * mag)

    fig = Figure(size = (1500, 520))

    ax1 = Axis(fig[1, 1], aspect = DataAspect(), title = "Analytical",   xlabel = "x", ylabel = "y")
    ax2 = Axis(fig[1, 3], aspect = DataAspect(), title = "Numerical",    xlabel = "x")
    ax3 = Axis(fig[1, 5], aspect = DataAspect(), title = "|Difference|", xlabel = "x")

    hm1 = heatmap!(ax1, gx, gy, an;   colormap = (Reverse(:matter), 1), crange...)
    hm2 = heatmap!(ax2, gx, gy, num;  colormap = (Reverse(:matter), 1), crange...)
    hm3 = heatmap!(ax3, gx, gy, diff; colormap = (Reverse(:matter), 1), colorrange = (-dmax, dmax))

    yx_ratio = (gy[end] - gy[1]) / (gx[end] - gx[1])
    rowsize!(fig.layout, 1, Aspect(1, yx_ratio))

    Colorbar(fig[1, 2], hm1, label = label, width = 14)
    Colorbar(fig[1, 4], hm2, label = label, width = 14)
    Colorbar(fig[1, 6], hm3, width = 14)

    hideydecorations!(ax2, grid = false)
    hideydecorations!(ax3, grid = false)
    colgap!(fig.layout, 8)

    for col in (2, 4, 6)
        colsize!(fig.layout, col, Fixed(70))
    end
    return fig
end

geometry_combo = GtkComboBoxText()
for g in ("elliptical", "circular")
    push!(geometry_combo, g)
end
geometry_combo.active = 0

preset_names = ("Custom", "Schmid2003", "Duretz2026_1_ps", "Duretz2026_2_ss",
                 "Duretz2026_3_exp", "Duretz2026_4_comp", "Duretz2026_5_mixed",
                 "Duretz2026_6_angled", "Duretz2026_7_oop")
preset_combo = GtkComboBoxText()
for p in preset_names
    push!(preset_combo, p)
end
preset_combo.active = length(preset_names) - 1   # Duretz2026_7_oop

field_names = ("P", "Vx", "Vy",
               "sigma_xx", "sigma_yy", "sigma_zz", "sigma_xy",
               "tau_xx", "tau_yy", "tau_zz", "tau_xy", "tau_II")
field_combo = GtkComboBoxText()
for f in field_names
    push!(field_combo, f)
end
field_combo.active = 0

entries = Dict{Symbol,Any}()
param_rows = (
    (:ηm,        "η_m  matrix viscosity",       "1.0"),
    (:ηi,        "η_i  inclusion viscosity",     "1e-2"),
    (:ξm,        "ξ_m  matrix bulk viscosity",   "1.0"),
    (:ξi,        "ξ_i  inclusion bulk viscosity","1.0"),
    (:ε̇dot,      "ε̇    pure shear rate",         "0.0"),
    (:γ̇dot,      "γ̇    simple shear rate",       "0.0"),
    (:ζ̇dot,      "ζ̇    out-of-plane shear rate", "1.0"),
    (:ε̇zzdot,    "ε̇zz  out-of-plane normal rate","1.0"),
    (:ri,        "ri   inclusion radius",        "0.1"),
    (:t,         "t    Joukowski shape param",   "2.0"),
    (:alpha_deg, "α    orientation (deg)",       "0.0"),
)
for (key, _, default) in param_rows
    entries[key] = GtkEntry(; text = default, hexpand = true)
end
nc_entry = GtkEntry(; text = "81", hexpand = true)

run_btn    = GtkButton("Run solver")
save_btn   = GtkButton("Save figure")
spinner    = GtkSpinner()
status_lbl = GtkLabel(""; wrap = true, hexpand = true, xalign = 0.0)

canvas = GtkCanvas(1500, 520; vexpand = true, hexpand = true)
current_fig = Ref{Any}(Figure())

@guarded draw(canvas) do widget
    screen = CairoMakie.Screen(current_fig[].scene, MAKIE_SCREEN_CONFIG, Gtk4.cairo_surface(canvas))
    CairoMakie.resize!(current_fig[].scene, Gtk4.width(widget), Gtk4.height(widget))
    CairoMakie.cairo_draw(screen, current_fig[].scene)
end

function refresh_plot!()
    state.ready || return
    current_fig[] = build_figure(active_text(field_combo))
    canvas.draw(canvas)
    reveal(canvas)
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
function getfloat(key::Symbol, default::Float64)
    v = tryparse(Float64, entries[key].text)
    if v === nothing
        entries[key].text = string(default)
        return default
    end
    return v
end

function apply_preset!(testname::String, geometry::Symbol)
    testname == "Custom" && return
    p = InclusionAnalytics.preset(Symbol(testname), geometry)
    entries[:ηm].text        = string(p.ηm)
    entries[:ηi].text        = string(p.ηi)
    entries[:ξm].text        = string(p.ξm)
    entries[:ξi].text        = string(p.ξi)
    entries[:ε̇dot].text      = string(p.ε̇)
    entries[:γ̇dot].text      = string(p.γ̇)
    entries[:ζ̇dot].text      = string(p.ζ̇)
    entries[:ε̇zzdot].text    = string(p.ε̇zz)
    entries[:ri].text        = "0.1"   # all current presets use ri = 0.1
    entries[:t].text         = string(p.t)
    entries[:alpha_deg].text = string(p.α * 180 / pi)
end

# ---------------------------------------------------------------------------
# Signals
# ---------------------------------------------------------------------------
function on_preset_or_geometry_changed(args...)
    @idle_add begin
        try
            apply_preset!(active_text(preset_combo), Symbol(active_text(geometry_combo)))
        catch e
            status_lbl.label = "Error applying preset: $(sprint(showerror, e))"
        end
        false
    end
    nothing
end

Gtk4.GLib.on_notify(on_preset_or_geometry_changed, preset_combo, "active")
Gtk4.GLib.on_notify(on_preset_or_geometry_changed, geometry_combo, "active")

Gtk4.GLib.on_notify(field_combo, "active") do args...
    @idle_add begin
        refresh_plot!()
        false
    end
    nothing
end

signal_connect(run_btn, "clicked") do widget
    geometry = Symbol(active_text(geometry_combo))
    ηm     = getfloat(:ηm, 1.0)
    ηi     = getfloat(:ηi, 1e-2)
    ξm     = getfloat(:ξm, 1.0)
    ξi     = getfloat(:ξi, 1.0)
    ε̇      = getfloat(:ε̇dot, 0.0)
    γ̇      = getfloat(:γ̇dot, 0.0)
    ζ̇      = getfloat(:ζ̇dot, 0.0)
    ε̇zz    = getfloat(:ε̇zzdot, 0.0)
    ri     = getfloat(:ri, 0.1)
    t      = getfloat(:t, 2.0)
    α      = getfloat(:alpha_deg, 0.0) * pi / 180
    nc_in  = tryparse(Int, nc_entry.text)
    nc     = nc_in === nothing ? 81 : clamp(nc_in, 11, 401)
    nc_entry.text = string(nc)

    run_btn.sensitive = false
    save_btn.sensitive = false
    start(spinner)
    status_lbl.label = "Running numerical solver (nc = $nc)... window will be unresponsive until done."
    flush_gtk_events!()   # make the message above actually paint before we block

    try
        params = build_params(geometry, ηm, ηi, ξm, ξi, ε̇, γ̇, ζ̇, ε̇zz, ri, t, α)
        f_anal = f_anal_for(geometry)
        V, P, σ, X    = InclusionAnalytics.Stokes2D(nc, :custom, params, f_anal)
        Va, Pa, σa    = InclusionAnalytics.analytics_stag(f_anal, X, params, geometry)
        Ve, Pe, σe    = InclusionAnalytics.errors_stag(V, P, σ, Va, Pa, σa)

        state.ready    = true
        state.geometry = geometry
        state.nc       = nc
        state.X, state.V, state.P, state.σ = X, V, P, σ
        state.Va, state.Pa, state.σa       = Va, Pa, σa
        state.Ve, state.Pe, state.σe       = Ve, Pe, σe
        refresh_plot!()
        status_lbl.label = @sprintf(
            "Done (nc=%d). mean|Vx err|=%.2e  mean|Vy err|=%.2e  mean|P err|=%.2e  mean|σxy err|=%.2e",
            nc, Statistics.mean(Ve.x), Statistics.mean(Ve.y), Statistics.mean(Pe), Statistics.mean(σe.xy))
    catch e
        status_lbl.label = "Error: $(sprint(showerror, e))"
    finally
        stop(spinner)
        run_btn.sensitive  = true
        save_btn.sensitive = true
    end
end

signal_connect(save_btn, "clicked") do widget
    state.ready || return
    outdir = joinpath(@__DIR__, "output")
    mkpath(outdir)
    fname = joinpath(outdir, "comparison_$(state.geometry)_nc$(state.nc)_$(active_text(field_combo)).png")
    CairoMakie.save(fname, current_fig[], px_per_unit = 2)
    status_lbl.label = "Saved: $fname"
end

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------
function build_param_grid()
    grid = GtkGrid(; column_spacing = 6, row_spacing = 6, margin_top = 8, margin_bottom = 8,
                     margin_start = 8, margin_end = 8)
    row = 1
    grid[1,row] = GtkLabel("Geometry"); grid[2,row] = geometry_combo; row += 1
    grid[1,row] = GtkLabel("Preset");   grid[2,row] = preset_combo;   row += 1
    for (key, label, _) in param_rows
        grid[1,row] = GtkLabel(label; xalign = 0.0)
        grid[2,row] = entries[key]
        row += 1
    end
    grid[1,row] = GtkLabel("nc  grid resolution"); grid[2,row] = nc_entry; row += 1
    grid[1,row] = run_btn; grid[2,row] = spinner; row += 1
    grid[1:2,row] = GtkLabel("Field to display"); row += 1
    grid[1:2,row] = field_combo; row += 1
    grid[1:2,row] = save_btn; row += 1
    grid[1:2,row] = status_lbl; row += 1
    return grid
end
grid = build_param_grid()

left_scroller = GtkScrolledWindow()
left_scroller.child = grid
left_scroller.min_content_width = 460   # wide enough for the longest label + entry, no horizontal scroll needed
left_scroller.hscrollbar_policy  = Gtk4.PolicyType_NEVER

right_box = GtkBox(:v)
push!(right_box, canvas)

main_box = GtkBox(:h; hexpand = true, vexpand = true)
push!(main_box, left_scroller)
push!(main_box, right_box)

win = GtkWindow(main_box, "InclusionAnalytics -- analytical vs numerical", 1900, 850)

apply_preset!(active_text(preset_combo), Symbol(active_text(geometry_combo)))

if !isinteractive()
    GLib.start_main_loop()
    GLib.waitforsignal(win, "close-request")
end

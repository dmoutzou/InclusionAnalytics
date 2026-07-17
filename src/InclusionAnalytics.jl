module InclusionAnalytics

using StaticArrays, Printf, Statistics, LinearAlgebra

include("Physical_parameters.jl")
export preset

include("AnalyticsCircularInclusion.jl")
export Analytics_circle

include("AnalyticsEllipticalInclusion.jl")
export joukowski, t_to_ri, ellipse_axes, to_zeta, Analytics_ellipse

include("Numerics.jl")
export Stokes2D, f_anal

end

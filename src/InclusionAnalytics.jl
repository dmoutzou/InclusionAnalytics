module InclusionAnalytics

using StaticArrays, Printf, Statistics, LinearAlgebra

include("Analysis.jl")
export preset_circle, preset_ellipse, preset

include("Numerics.jl")
export Stokes2D

include("AnalyticsCircularInclusion.jl")
export Analytics_old, Analytics_new

include("AnalyticsEllipticalInclusion.jl")
export joukowski, t_to_ri

end # module InclusionAnalytics

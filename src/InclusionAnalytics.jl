module InclusionAnalytics

using StaticArrays, Printf, Statistics, LinearAlgebra

include("PhysicalParameters.jl")
export preset

include("AnalyticsCircularInclusion.jl")
export analytics_circle

include("AnalyticsEllipticalInclusion.jl")
export joukowski, t_to_ri, ellipse_axes, to_zeta, analytics_ellipse

include("Numerics.jl")
export Stokes2D, f_anal

include("AnalyticsStagerredGrid.jl")
export analytics_stag, errors_stag

end

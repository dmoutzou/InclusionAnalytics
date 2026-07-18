module InclusionAnalytics

using StaticArrays, Printf, Statistics, LinearAlgebra

include("PhysicalParameters.jl")
export preset

include("Analytics.jl")
export f_anal, analytics_circle, analytics_ellipse

include("Numerics.jl")
export Stokes2D, f_anal

include("AnalyticsStagerredGrid.jl")
export analytics_stag, errors_stag

end

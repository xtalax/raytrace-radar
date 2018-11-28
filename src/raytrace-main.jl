include("raytracer-types.jl")
include("Surface.jl")

function Fresnel(ray::Ray, interaction::Intersection)
    

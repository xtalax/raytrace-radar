module TriangleIntersect

export Point, Triangle, Ray, intersect, Intersection, no_intersection

const no_intersection = Intersection(Point(0,0,0),Point(0,0,0) 0.0, false)

function triangleIntersect(r::Ray, t::Triangle)
    denom = t.normal*r.direction
    denom == 0 && return no_intersection
    ri = t.normal*(t. - r.origin) / denom
    ri <= 0 && return no_intersection
    plane_intersection =  ri * r.direction + r.origin
    w = plane_intersection - t.a
    wv1 = w*t.v1
    wv2 = w*t.v2
    s_intersection = (t.v1v2*wv2 - t.v2v2*wv1) / t.denom
    s_intersection <= 0 && return no_intersection
    s_intersection >= 1 && return no_intersection
    t_intersection = (t.v1v2*wv1 - t.v1v1*wv2) / t.denom
    t_intersection <= 0 && return no_intersection
    t_intersection >= 1 && return no_intersection
    s_intersection + t_intersection >= 1 && return no_intersection
    return Intersection(t.a + s_intersection*t.v1+t_intersection*t.v2,t.normal, ri, true)
end

end # module

const no_intersection = Intersection(Point(0,0,0), 0.0, Point(0,0,0), false)

function polygonIntersect(r::Ray, t::Triangle)
    denom = t.normal*r.direction
    denom == 0 && return no_intersection
    ri = t.normal*(t.A.Coords - r.origin) / denom
    ri <= 0 && return no_intersection
    plane_intersection =  ri * r.direction + r.origin
    w = plane_intersection - t.A.Coords
    wv1 = w*t.v1
    wv2 = w*t.v2
    s_intersection = (t.v1v2*wv2 - t.v2v2*wv1) / t.denom
    s_intersection <= 0 && return no_intersection
    s_intersection >= 1 && return no_intersection
    t_intersection = (t.v1v2*wv1 - t.v1v1*wv2) / t.denom
    t_intersection <= 0 && return no_intersection
    t_intersection >= 1 && return no_intersection
    s_intersection + t_intersection >= 1 && return no_intersection

    return Intersection(t.A.Coords + s_intersection*t.v1+t_intersection*t.v2, ri,t.normal, true)
end

function polygonIntersect(r::Ray, p::Parallelogram)
    denom = p.normal*r.direction
    denom == 0 && return no_intersection
    ri = p.normal*(p.origin - r.origin) / denom
    ri <= 0 && return no_intersection
    plane_intersection =  ri * r.direction + r.origin
    w = plane_intersection - p.origin
    q1 = project(p.v1, w) ; q2 = project(p.v2, w)
    0 ≤ norm(q1) ≤ norm(p.v1) || return no_intersection
    0 ≤ norm(q2) ≤ norm(p.v2) || return no_intersection
    
    return Intersection(plane_intersection, ri, p.normal, true)
end

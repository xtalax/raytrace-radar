

function setreduce(x::Vector{T}, y::Vector{T}) where T <: Number
    x_out::Vector{T} = zeros(1)
    y_out::Vector{T} = zeros(1)

    for i in eachindex(x)
        if x[i] ∈ x_out
            y_out[x_out .== x[i]] .+= y[i]
        else
            push!(x_out, x[i])
            push!(y_out, y[i])
        end
    end

    return x_out, y_out
end
function duplicateflatten(detector)
    timess =detector.pixels[1].times
    timesp =detector.pixels[1].times
    amps= detector.pixels[1].ampS
    ampp =detector.pixels[1].ampP

    times, ampS = setreduce(timess, amps)
    times, ampP = setreduce(timesp, ampp)
    return times, ampS, ampP
end

function flattenarray(arr)
    rst = Any[]
    grep(v) = for x in v
        if isa(x, Array) grep(x) else push!(rst, x) end
    end
    grep(arr)
    rst
end

function nthparent(object, parentname, n=1::Int)
    return eval(Symbol(:object, String(:parentname)^n))
end

function raytraverse(ray::Ray)
    interaction_points::Vector{Point} = []
    while true
        if typeof(nthparent(ray, parent, i)) <: Ray
            push!(interaction_points, nthparent(ray, parent, i).origin)
        else
            return interaction_points
        end
    end
end

function linecreate(interaction_points::Vector{Point})

    lines::Vector{Vector{Point}} = [[]]
    length(interaction_points) = npoints
    sizehint!(lines, npoints - 1)

    for p in 1:npoints-1
        push!(lines, [interaction_points[p], interaction_points[p+1]])
    end
    return lines
end

function plotprepare(points::Vector{Point})
    for p in points
        push!(x, p.x)
        pusy!(y, p.y)
        push!(z, p.z)
    end
    cartesian = ["x" => x,
                 "y" => y,
                 "z" => z]
    return cartesian
end

function foo(origin::Point, x::Point, y::Point)
    shape = Parallelogram(origin,x,y)
    times::Vector{Float64} = [0.0]
    ampS::Vector{Complex} = []
    ampP::Vector{Complex}= []
    ϕs::Vector{Complex}= []
    ϕp::Vector{Complex} =[]
    return shape, times, ampS, ampP, ϕs, ϕp
end

convert

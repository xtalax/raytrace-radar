using LinearAlgebra
using StaticArrays
using Unitful #uconvert() to convert units
using Rotations
#using PhysicalConstants


const c₀ = 299792458u"m/s"
const ε₀ = (8.85418782*10^-12)u"F/m"
const μ₀ = (1.2566370614*10^-6)u"H/m"

## Utils
unitize(Vec::AbstractVector{Number}) = Vec/norm(Vec)


rad(ϴ::typeof(0.0u"rad")) = ϴ
rad(ϴ::typeof(0.0u"°")) = uconvert(u"rad",ϴ)
rad(ϴ::Float64) = (ϴ)u"rad"
°(ϴ::typeof(0.0u"°")) = ϴ
°(ϴ::typeof(0.0u"rad")) = uconvert(u"°",ϴ)
°(ϴ::Real) = (ϴ)u"°"


abstract type Geometry end
abstract type Shape <: Geometry end
abstract type Polygon <: Shape end

# Point3D type

mutable struct Point<:Geometry
    x::Union{Float64,typeof(1.0u"m")}
    y::Union{Float64,typeof(1.0u"m")}
    z::Union{Float64,typeof(1.0u"m")}
    Point(x::Float64, y::Float64, z::Float64) = new(x,y,z)
    Point(A::Union{Vector{Float64},SVector{3,Float64}})=length(A)==3 ?  new(A[1],A[2],A[3]) : error("Point must be instantiated with 3 floats or an array of 3 floats")

end
import Base.*
import Base./
import Base.+
import Base.-
import LinearAlgebra.cross

*(p::Point, n::Number) = Point(p.x*n, p.y*n, p.z*n)
*(n::Number, p::Point) = p*n
*(p1::Point, p2::Point) = p1.x*p2.x+p1.y*p2.y+p1.z*p2.z
-(p1::Point, p2::Point) = Point(p1.x-p2.x, p1.y-p2.y, p1.z-p2.z)
+(p1::Point, p2::Point) = Point(p1.x+p2.x, p1.y+p2.y, p1.z+p2.z)
cross(p1::Point, p2::Point) = Point(p1.y*p2.z-p1.z*p2.y, -p1.x*p2.z+p1.z*p2.x, p1.x*p2.y-p1.y*p2.x)
/(p::Point, n::Number) = Point(p.x/n, p.y/n, p.z/n)
unitize(p::Point) = p/(p*p)

#### Defines the node class - Use the buildNode method to construct and track already existing nodes properly



mutable struct Node <:Geometry
    Coords::Point

    InEdges::Array{Geometry,1} #all the lines that this node defines
    InTris::Array{Polygon,1}    #all the triangles that this node defines

end
function buildNode(Coords,NodeList)
    ExtantNode = get(NodeList,Coords,false)
    if ExtantNode == false
        self = Node(Coords,[],[])
        push!(NodeList,Coords => self)
        return self
    end
    return ExtantNode
end
#### Defines the Edge class - Use the buildEdge method to construct and track already existing Edges properly

mutable struct Edge <:Geometry
    A::Node
    B::Node
    InTris::Array{Polygon,1}
    # needs a constructor to handle only having Point3Ds passed to it
    Edge(A,B)=new(A,B,[])
end

function buildEdge(A::Node, B::Node,EdgeList)

    ExtantEdge = get(EdgeList,[A.Coords,B.Coords],false)
    Permuted = get(EdgeList,[B.Coords,A.Coords],false)
    if ExtantEdge == false
        if Permuted == false
            self = Edge(A,B); push!(EdgeList,[A.Coords,B.Coords] => self)
            push!(A.InEdges, self) ; push!(B.InEdges, self)
            return self
        end
        return Permuted
    end
    return ExtantEdge

end

#### Defines the Triangle class - Use the buildEdge method to construct and track already existing Edges properly
######## This class turned out to be a monster
mutable struct Triangle <: Polygon

    A::Node
    B::Node
    C::Node
    AB::Edge
    BC::Edge
    CA::Edge
    v1::Point
    v2::Point
    normal::Point
    v1v1::Float64
    v2v2::Float64
    v1v2::Float64
    denom::Float64
    center::Point
    area::Float64

end

function buildTriangle(A::Node,B::Node,C::Node,AB::Edge,BC::Edge,CA::Edge)
    a = A.Coords
    b = B.Coords
    c = C.Coords

    v1 = b-a
    v2 = c-a
    normal = cross(v1, v2)
    v1v1 = v1*v1
    v2v2 = v2*v2
    v1v2 = v1*v2
    denom = v1v2*v1v2 - v1v1*v2v2
    center = (a+b+c)/3
    area = normal/2
    self = Triangle(A,B,C,AB,BC,CA, v1, v2, normal, v1v1, v2v2, v1v2, denom,center,area)
    push!(A.InTris, self) ; push!(B.InTris, self) ;  push!(C.InTris, self)
    push!(AB.InTris, self) ; push!(BC.InTris, self) ;  push!(CA.InTris, self)
    return self

end
buildTriangle(Points::Array{Point,1}) = buildTriangle(Points[1], Points[2], Points[3])
buildTriangle(Nodes::Array{Node,1}) = buildTriangle(Nodes[1], Node[2], Node[3])





struct Surface <: Geometry
    Triangles::Array{Triangle,1}
    Edges::Array{Edge,1}
    Nodes::Array{Node,1}
    centroid::Point
    extent::Float64
    RefrIndex::Float64
    function Surface(nodeCoords,triIndices,nTriangles,RefrIndex,origin = [0.0,0.0,0.0])
        Triangles::Array{Triangle,1} = []
        Edges::Array{Edge,1} = []
        Nodes::Array{Node,1}= []
        NodeList::Dict = Dict([])
        EdgeList::Dict = Dict([])
        extent::Float64
        for i in nTriangles

            A = buildNode(Point(nodeCoords[:,triIndices[1,i]]),NodeList)
            B = buildNode(Point(nodeCoords[:,triIndices[2,i]]),NodeList)
            C = buildNode(Point(nodeCoords[:,triIndices[3,i]]),NodeList)
            ra = sqrt(sum(nodeCoords[:,triIndices[1,i]]).^2))
            rb = sqrt(sum(nodeCoords[:,triIndices[2,i]]).^2))
            rc = sqrt(sum(nodeCoords[:,triIndices[3,i]]).^2))
            extent = max([ra,rb,rc,extent])
            pointsum = nodeCoords[:,triIndices[1,i]])+nodeCoords[:,triIndices[2,i]])+nodeCoords[:,triIndices[3,i]])

            AB = buildEdge(A,B,EdgeList)
            BC = buildEdge(B,C,EdgeList)
            CA = buildEdge(C,A,EdgeList)
            ABC = buildTriangle(A, B, C, AB, BC, CA)
            push!(Triangles,ABC)
            push!(Edges, AB,BC,CA)
            push!(Nodes, A, B, C)
        end
        centroid = pointsum./(nTriangles*3)
        return new(Triangles,Edges,Nodes,centroid,extent,RefrIndex)
    end
end

mutable struct Scene
    

struct Intersection
    point::Point # intersecting point
    distance::Float64 # intersecting distance
    normal::Point
    isHit::Bool
end

abstract type Beam

mutable struct Ray <: Beam#The Rays object
    origin::Vector{typeof(1.0u"m")}
    direction::Vector{Real}

    parent::Beam
    decendents::Vector{Beam}
    isActive::Bool
    tof::typeof(1.0u"s") #time of flight of the Ray
    intensity::typeof(1.0u"J")
    c::typeof(1.0u"m/s")
    function Ray(origin, direction, parent tof=0.0u"s", S=1.0u"J", c=c₀)


          new(origin, direction,parent,[]::Vector{Beam},true tof, S, c)

     end

end

struct Source <: Beam
    ray::Ray
    intensity::typeof(1.0u"J")
    isActive::Bool
end
#=
Bundle returns an array of rays. it takes a coord system origin, an direction,
an angle to rotate about that direction, the number of rays along a side,
the shape of the generated bundle, and 2 length parameters
shape = "square"


=#
mutable struct Bundle
    Rays::Vector{Ray}
    N::Int64
    function Bundle(origin = repeat([0.0u"m"], 3)::Vector{typeof(1.0u"m")}, direction=[0.0,0.0,1.0]::Vector{Real}, φ = 0, N = 1000::Int, shape = "rect"::String, X₁=1.0u"m"::typeof(u"m"), X₂=1.0u"m"::typeof(u"m"))
        if shape == "rect"
            A = X₂/X₁
            r_z=RotZ(rad(φ)/u"rad") # rotates about the z axis by angle φ
            q=rotation_between([0.0,0.0,1.0],normalize!(direction)) #generates a quaternion matrix that rotates Point3Ds an angle and direction between the args about the origin
            Origins=[Matrix(q*r_z)*[X₁*x, X₂*y, 0.0u"m"]-origin   for x in -0.5:1*A/(sqrt(N)-1):0.5, y in -0.5:1/(sqrt(N)*A-1):0.5] #a hench comprehension that generates a square of evenly spaced origins
            N = size(Origins,1)*size(Origins,2)

            Origins = reshape(Origins, (N, ))
            directions = [direction for ii in 1:N]
            rays = Ray.(Origins, directions)
            new(rays,N)
        elseif shape == "sphere"
            new(rays,N)
        end
    end
end






r = Ray([1.0u"m", 2.0u"m", 3.0u"m"],[1,0,0])


B=Bundle([1.0, 2.0, 3.0]u"m", [0.0, -1.0, -1.0], π/2, 500, 1.0u"m")

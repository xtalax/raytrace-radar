

#=export
Geometry, Shape, Polygon, Point, Node, Edge, Triangle, Surface, Parallelogram, Intersection, Beam, Source, Detector, Ray, Bundle,
LinearAlgebra.⋅, Base.*, LinearAlgebra.×, Base.+, Base.-, LinearAlgebra.Vector,  Base.isnan, LinearAlgebra.normalize, LinearAlgebra.norm
end=#

const c₀ = 299792458
const ε₀ = (8.85418782*10^-12)
const μ₀ = (1.2566370614*10^-6)

## Utils
unitize(Vec::AbstractVector{Number}) = Vec/norm(Vec)





abstract type Geometry end
abstract type Shape <: Geometry end
abstract type Polygon <: Shape end

# Point3D type

struct Point<:Geometry
    x::Float64
    y::Float64
    z::Float64
    Point(A::Array{Float64,1})=length(A)==3 ?  new(A[1],A[2],A[3]) : error("Arrays converted to Points must be of length 3")
    Point(A::SVector{3,Float64}) = new(A[1],A[2],A[3])
    Point(x::Int64, y::Int64, z::Int64) =  new(Float64(x),Float64(y),Float64(z))
    Point(x::Number,y::Number,z::Number) = new(Float64(x),Float64(y),Float64(z))
end

import Base.*
import Base./
import Base.+
import Base.-
import Base.^
import Base.isnan
import LinearAlgebra.cross
import LinearAlgebra.normalize
import LinearAlgebra.⋅
import LinearAlgebra.norm


*(p::Point, n::Number) = Point(p.x*n, p.y*n, p.z*n)
*(n::Number, p::Point) = p*n
*(p1::Point, p2::Point) = p1.x*p2.x+p1.y*p2.y+p1.z*p2.z
⋅(a::Point,b::Number) = *(a,b)
⋅(a::Number,b::Point) = *(a,b)
⋅(a::Point,b::Point) = *(a,b)
-(p1::Point, p2::Point) = Point(p1.x-p2.x, p1.y-p2.y, p1.z-p2.z)
-(p::Point) = Point(-p.x,-p.y,-p.z)
+(p1::Point, p2::Point) = Point(p1.x+p2.x, p1.y+p2.y, p1.z+p2.z)
Vector(p::Point) = [p.x,p.y,p.z]
SVector(p::Point) = SVector(p.x,p.y,p.z)
cross(p1::Point, p2::Point) = Point(p1.y*p2.z-p1.z*p2.y, -p1.x*p2.z+p1.z*p2.x, p1.x*p2.y-p1.y*p2.x)
isnan(p::Point) = isnan(p.x) | isnan(p.y) | isnan(p.z)
/(p::Point, n::Number) = Point(p.x/n, p.y/n, p.z/n)
norm(p::Point) = sqrt(p.x^2+p.y^2+p.z^2)

unitize(p::Point) = p/norm(p)
normalize(p::Point) = unitize(p)
normalize(c::Complex{Float64}) =  c/norm(c)


project(b::Point, a::Point) = (a*unitize(b))*unitize(b)

vectorAngle(p1::Point, p2::Point) = acos(p1*p2/(norm(p1)*norm(p2)))
vectorAngle(p1::SVector{3,Float64}, p2::SVector{3,Float64}) = SVector(vectorAngle(Point(p1),Point(p2)))
^(p::Point,n::Number) = Point(p.x^n,p.y^n,p.z^n)
#### Defines the node class - Use the buildNode method to construct and track already existing nodes properly



mutable struct Node <:Geometry
    Coords::Point

    InEdges::Array{Geometry,1} #all the lines that this node defines
    InTris::Array{Polygon,1}    #all the triangles that this node defines

end
function buildNode(Coords,NodeList = Dict([]))
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

function buildEdge(A::Node, B::Node,EdgeList=Dict([]))

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

    normal = v1×v2
    v1v1 = v1*v1
    v2v2 = v2*v2
    v1v2 = v1*v2
    denom = v1v2*v1v2 - v1v1*v2v2
    center = (a+b+c)/3
    area = norm(normal)/2
    normal = normalize(normal)
    self = Triangle(A,B,C,AB,BC,CA, v1, v2, normal, v1v1, v2v2, v1v2, denom,center,area)

    push!(A.InTris, self) ; push!(B.InTris, self) ;  push!(C.InTris, self)
    push!(AB.InTris, self) ; push!(BC.InTris, self) ;  push!(CA.InTris, self)
    return self

end
function buildTriangle(Points::Array{Point,1})
    A,B,C = buildNode(Points[1]), buildNode(Points[2]), buildNode(Points[3])
    AB,BC,CA = buildEdge(A,B),buildEdge(B,C),buildEdge(C,A)
    Triangle(A,B,C,AB,BC,CA)
end

buildTriangle(Nodes::Array{Node,1}) = buildTriangle(Nodes[1], Node[2], Node[3])

struct Parallelogram <:Polygon
    origin::Point
    v1::Point
    v2::Point
    normal::Point
    function Parallelogram(origin::Point,v1::Point,v2::Point)
        normal = normalize(cross(v1,v2))

        return new(origin, v1, v2, normal)
    end
end

mutable struct Surface <: Geometry
    Triangles::Array{Triangle,1}
    Edges::Array{Edge,1}
    Nodes::Array{Node,1}
    centroid::Point
    extent::Float64
    RefrIndex::Float64
    function Surface(nodeCoords,triIndices,nTriangles, RefrIndex=1.0::Float64, scale=1.0::Float64,origin = [0.0,0.0,0.0])
        Triangles::Array{Triangle,1} = []
        Edges::Array{Edge,1} = []
        Nodes::Array{Node,1}= []
        NodeList::Dict = Dict([])
        EdgeList::Dict = Dict([])
        extent::Float64 = 0.0
        pointsum = SVector(0.0,0.0,0.0)*0.0
        for i in 1:nTriangles

            A = buildNode(Point(nodeCoords[:,triIndices[1,i]]*scale),NodeList)
            B = buildNode(Point(nodeCoords[:,triIndices[2,i]]*scale),NodeList)
            C = buildNode(Point(nodeCoords[:,triIndices[3,i]]*scale),NodeList)

            pointsum::SVector{3,Float64} = (nodeCoords[:,triIndices[1,i]]+nodeCoords[:,triIndices[2,i]]+nodeCoords[:,triIndices[3,i]])*scale + pointsum

            AB = buildEdge(A,B,EdgeList)
            BC = buildEdge(B,C,EdgeList)
            CA = buildEdge(C,A,EdgeList)

            ABC = buildTriangle(A, B, C, AB, BC, CA)

            push!(Triangles,ABC)
            push!(Edges, AB,BC,CA)
            push!(Nodes, A, B, C)

        end
        centroid = Point(pointsum./(nTriangles*3))
        extent = maximum(norm(node.Coords-centroid) for node in Nodes)
        new(Triangles,Edges,Nodes,centroid,extent,RefrIndex)
    end
end

mutable struct Pixel <:Geometry
    shape::Parallelogram
    TEBasis::Point
    times::Vector{Float64}
    ampS::Vector{Float64}
    ampP::Vector{Float64}
    ϕs::Vector{Complex}
    ϕp::Vector{Complex}

    function Pixel(origin::Point,x::Point,y::Point,TEBasis::Point)
        shape = Parallelogram(origin,x,y)
        times::Vector{Float64} = []
        ampS::Vector{Complex} = []
        ampP::Vector{Complex}= []
        ϕs::Vector{Complex}= []
        ϕp::Vector{Complex} =[]
        return new(shape,TEBasis, times, ampS, ampP,ϕs, ϕp)
    end
end


mutable struct Detector <:Geometry
    pixels::Vector{Pixel}
    shape::Polygon
    function Detector(origin = Point(0.0,0.0,0.0)::Point, direction=Point([0.0,0.0,1.0])::Point, φ = 0.0, Nx=1::Int, Ny=1::Int, shape = "square"::String, X=1.0, Y=1.0)

        Yvec = [0.0, Y, 0.0]
        Xvec = [X, 0.0, 0.0]
        TEvec =[0.0,1.0,0.0]

        r_z=RotZ(φ) # rotates about the z axis by angle φ
        q=rotation_between([0.0,0.0,1.0],Vector(normalize(direction))) #generates a quaternion matrix that rotates Point3Ds an angle and direction between the args about the origin
        T = Matrix(q*r_z)

        Ypoint = Point(T*Yvec)
        Xpoint = Point(T*Xvec)
        TEBasis = Point(T*TEvec)

        y = Ypoint/Ny
        x = Xpoint/Nx

        shapeOrigin = origin-Ypoint/2-Xpoint/2
        Origins=[shapeOrigin+x*nx+y*ny   for nx in 0:Nx-1, ny in 0:Ny-1]
        N = Nx*Ny

        Origins = reshape(Origins, (N, ))

        pixels = Pixel.(Origins, (x for i in 1:N), (y for i in 1:N), (TEBasis for i in 1:N))
        shape = Parallelogram(shapeOrigin, Xpoint, Ypoint)
        new(pixels,shape)

    end
end
struct Intersection
    point::Point # intersecting point
    distance::Float64 # intersecting distance
    normal::Point
    isHit::Bool
end

abstract type Beam end

struct Ray <: Beam#The Rays object
    budget::Int64
    origin::Point
    direction::Point
    TEBasis::Point
    parent::Beam
    ξ::Float64 #cumulative optical path length of the ray
    ampS::Float64
    ampP::Float64
    ϕs::Complex{Float64}
    ϕp::Complex{Float64}
    n::Float64
    function Ray(budget, origin, direction, TEBasis, parent, ξ=0.0, ampS = 1/sqrt(2), ampP = 1/sqrt(2), ϕs =1.0+0.0im, ϕp=1.0+0.0im,  n=1)


          new(budget, origin, normalize(direction), normalize(TEBasis) ,parent, ξ, ampS, ampP, ϕs, ϕp,  n)

     end

end

struct Source <: Beam

    intensity::Float64
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
    function Bundle(origin = Point(0.0,0.0,0.0)::Point, direction=Point([0.0,0.0,1.0])::Point, φ = 0.0, Nx=60::Int, Ny=60::Int, shape = "square"::String, X=0.6, Y=0.6)

        Yvec = [0.0, Y, 0.0]
        Xvec = [X, 0.0, 0.0]
        TEBasis =[0.0,1.0,0.0]

        r_z=RotZ(φ) # rotates about the z axis by angle φ
        q=rotation_between([0.0,0.0,1.0],Vector(normalize(direction))) #generates a quaternion matrix that rotates Point3Ds an angle and direction between the args about the origin
        T = Matrix(q*r_z)
        Ypoint = Point(T*Yvec)
        Xpoint = Point(T*Xvec)
        TEBasis = Point(T*TEBasis)
        y = Ypoint/Ny
        x = Xpoint/Nx

        shapeOrigin = origin-Ypoint/2-Xpoint/2
        Origins=[shapeOrigin+x*nx+y*ny + y/2 + x/2   for nx in 0:Nx-1, ny in 0:Ny-1]
        N = Nx*Ny

        Origins = reshape(Origins, (N, ))
        rays = Ray.(20, Origins, (direction for i in 1:N), (TEBasis for i in 1:N), (Source(1.0,true) for i in 1:N))

        new(rays, N)

    end
end

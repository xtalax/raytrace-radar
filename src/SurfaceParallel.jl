
const nullIntersection = Intersection(Point(0,0,0), 0.0, Point(0,0,0), false)
@everywhere const GPUnullIntersection = [0.0,0.0,0.0, 0.0, 0.0,0.0,0.0, 0]
const dummyray = Ray(Point(0.0,0.0,0.0),Point(0.0,1.0,0.0), Point(1.0,0.0,0.0), Source(1.0, true))

function polygonIntersect(r::Ray, t::Triangle)
    normal = r.direction⋅t.normal < 0 ? t.normal : -t.normal

    denom = normal*r.direction
    denom == 0 && return nullIntersection
    ri = normal*(t.A.Coords - r.origin) / denom
    ri <= 0 && return nullIntersection

    plane_intersection =  ri * r.direction + r.origin
    w = plane_intersection - t.A.Coords
    wv1 = w*t.v1
    wv2 = w*t.v2
    s_intersection = (t.v1v2*wv2 - t.v2v2*wv1) / t.denom
    s_intersection <= 0 && return nullIntersection
    s_intersection >= 1 && return nullIntersection
    t_intersection = (t.v1v2*wv1 - t.v1v1*wv2) / t.denom
    t_intersection <= 0 && return nullIntersection
    t_intersection >= 1 && return nullIntersection
    s_intersection + t_intersection >= 1 && return nullIntersection

    P = t.A.Coords + s_intersection*t.v1+t_intersection*t.v2
    r = norm(P-r.origin)
    return Intersection(P, ri, normal, true)
end



@everywhere function GPUIntersect(A::Vector{Float64}, normal::Vector{Float64}, v1::Vector{Float64}, v2::Vector{Float64}, v1v1::Float64, v2v2::Float64, v1v2::Float64, denom::Float64, direction::Vector{Float64}, origin::Vector{Float64})
    #throw("we made it reddit")
    denoma = normal⋅direction
    denoma == 0 && return GPUnullIntersection
    ri = normal⋅(A .- origin) ./ denoma
    ri <= 0 && return GPUnullIntersection

    plane_intersection =  ri .* direction .+ origin
    w = plane_intersection .- A
    wv1 = w⋅v1
    wv2 = w⋅v2
    s_intersection = (v1v2*wv2 - v2v2*wv1) / denom
    s_intersection <= 0 && return GPUnullIntersection
    s_intersection >= 1 && return GPUnullIntersection
    t_intersection = (v1v2*wv1 - v1v1*wv2) / denom
    t_intersection <= 0 && return GPUnullIntersection
    t_intersection >= 1 && return GPUnullIntersection
    s_intersection + t_intersection >= 1 && return GPUnullIntersection

    P = A .+ s_intersection.*v1.+t_intersection.*v2
    @show P
    r = norm(P.-origin)
    return [P[1],P[2],P[3], ri, normal[1], normal[2], normal[3], 1.0]
end

function polygonIntersect(r::Ray, p::Parallelogram)
    normal = r.direction⋅p.normal < 0 ? p.normal : -p.normal

    denom = normal*r.direction
    denom == 0 && return nullIntersection
    ri = normal*(p.origin - r.origin) / denom
    ri <= 0 && return nullIntersection
    plane_intersection =  ri * r.direction + r.origin
    w = plane_intersection - p.origin
    q1 = project(p.v1, w) ; q2 = project(p.v2, w)
    0 ≤ norm(q1) ≤ norm(p.v1) || return nullIntersection
    0 ≤ norm(q2) ≤ norm(p.v2) || return nullIntersection
    return Intersection(plane_intersection, ri, normal, true)
end



#### Function to allocate correct range to workers
@everywhere function getrange(n)
    @show tid = myid()-1
    nt = nprocs()-1
    d , r = divrem(n, nt)
    from = (tid - 1) * d + min(r, tid - 1) + 1
    to = from + d - 1 + (tid ≤ r ? 1 : 0)
    @show from:to
    from:to
end

#=
function run_bench(A, normal, v1v1, v2v2, v1v2, denom , direction, origin, out)

  # use dot syntax to apply `juliaset` to each elemt of q_converted
  # and write the output to result
  out .= polygonIntersect.(A, normal, v1v1, v2v2, v1v2, denom , direction, origin, ray)
  # all calls to the GPU are scheduled asynchronous,
  # so we need to synchronize
  GPUArrays.synchronize(out)
end=#

function SurfaceIntersect(surface::Surface, ray::Ray) # Outputs the index of the intersected triangle and the point of intersection.
    #(ray.direction*(ray.origin-surface.centroid))^2-norm(ray.origin-surface.centroid)^2-surface.extent^2>0 && return  nullIntersection,0 #Check if the ray goes anywhere near the sphere
    nIntersect = 0
    closestIntersection = nullIntersection


    #if length(surface.Triangles) < 0 # change later for gpu implementation
 ##### unpack tris and rays ########
      A = [Vector(t.A.Coords) for t in surface.Triangles]
      normal =[ Vector(ray.direction⋅t.normal < 0 ? t.normal : -t.normal) for t in surface.Triangles]
      v1 = [Vector(t.v1) for t in surface.Triangles]
      v2 = [Vector(t.v2) for t in surface.Triangles]
      v1v1 = [t.v1v2 for t in surface.Triangles]
      v2v2 = [t.v2v2 for t in surface.Triangles]
      v1v2= [t.v1v2 for t in surface.Triangles]
      denom = [t.denom for t in surface.Triangles]
      origin = Vector(ray.origin)
      direction = Vector(ray.direction)
      sendto(workers(), A=A)
      sendto(workers(), normal = normal)
      sendto(workers(), v1 = v1)
      sendto(workers(), v2 = v2)
      sendto(workers(), v1v1=v1v1)
      sendto(workers(), v2v2 = v2v2)
      sendto(workers(), v1v2 = v1v2)
      sendto(workers(), denom=denom)
      sendto(workers(), origin=origin)
      sendto(workers(), direction=direction)


#

#    else

      nT = length(triangles)

      intersections= SharedArray(Float64.(zeros(8,nT)))
      @passobj 1 workers() intersections


      println(length(A),length(normal),length(v1),length(v2),length(v1v2), length(v2v2),length(v1v2),length(denom), " ", size(intersections))
      @time @sync @distributed for k in 1:nprocs()-1
          range = getrange(nT)
           for i in range
               @show i
               intersections[:,i] .= GPUIntersect(A[i], normal[i], v1[i], v2[i], v1v1[i], v2v2[i], v1v2[i], denom[i], direction, origin)


           end

       end
      #=@time @sync @distributed for k in 1:nprocs()-1
          throw("im in ur loop throwin ur errors")
          range = getrange(nT)
          for i in range
              #@show A[i]
              #@show normal[i]
              #@show v1[i]
              #@show v2[i]
              #@show v1v2[i]
              #@show v2v2[i]
              #@show v1v2[i]
              #ntersections[:,i] .= GPUIntersect(A[i], normal[i], v1[i], v2[i], v1v1[i], v2v2[i], v1v2[i], denom[i], direction, origin)
              throw("winner winner")
          end

      end=#
      a = Array(intersections)
      I_D = a[argmin(a[4,:])]
      Intersection(Point(I_D[1],I_D[2],I_D[3]), I_D[4], Point(I_D[5],I_D[6],I_D[7]), Bool(I_D[8]))
      nIntersect =Int64(sum(I_D[8] for i in a))
    #  end

  return closestIntersection, nIntersect
end




@everywhere function load_ply_file(fileName::String, RefrIndex=1.0::Float64, scale=1.0::Float64, origin = [0.0,0.0,0.0]::Vector{Float64})
  println(" ~~~ loading surface mesh ~~~")
  nNodes = 0
  nTriangles = 0
  iHeader = 0
  i = 0
  iFile = open(fileName, "r")
  while !eof(iFile)
    line = readline(iFile)
    if occursin( "element vertex",line)
      nNodes = parse(Int, split(line, " ")[3])
    elseif occursin("element face",line)
      nTriangles = parse(Int, split(line, " ")[3])
    elseif occursin( "end_header",line)
      iHeader = i
      break
    end
    i += 1
  end
  close(iFile)

  println("       nNodes     : ", nNodes)
  println("       nTriangles : ", nTriangles)
  println("       iHeader    : ", iHeader)


  triIndices = zeros(Int64, 3, nTriangles)
  nodeCoords = zeros(Float64,3,nNodes)
  i = 0
  iFile = open(fileName, "r")
  while !eof(iFile)
    line = readline(iFile)
    if iHeader < i <= iHeader+nNodes

      xyz=collect(m.match for m in eachmatch(r"-?\d+(\.\d+)?", line))

      nodeCoords[1,i-iHeader] = parse(Float64,xyz[1])
      nodeCoords[2,i-iHeader]= parse(Float64,xyz[2])
      nodeCoords[3,i-iHeader]= parse(Float64,xyz[3])

    elseif i > iHeader+nNodes
      ijk=xyz=collect(m.match for m in eachmatch(r"(\d+)", line))

      triIndices[1,i-iHeader-nNodes] = parse(Int, ijk[2])+1
      triIndices[2,i-iHeader-nNodes]= parse(Int, ijk[3])+1
      triIndices[3,i-iHeader-nNodes] = parse(Int, ijk[4])+1
    end
    i += 1
  end
  close(iFile)
  println("~~~ Done ~~~ \n\n~~~ Building Surface ~~~")
  println(nodeCoords)
  println("tri indicies")
  println(triIndices)
  println(nTriangles)
  surface = Surface(nodeCoords, triIndices, nTriangles, RefrIndex)



  totalSurfaceArea = sum([triangle.area for triangle in surface.Triangles])
  return surface
end

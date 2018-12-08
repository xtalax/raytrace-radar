
const nullIntersection = Intersection(Point(0,0,0), 0.0, Point(0,0,0), false)

function SurfaceIntersect(surface::Surface, ray::Ray) # Outputs the index of the intersected triangle and the point of intersection.
    lIntersectMin = 0.0
    #(ray.direction*(ray.origin-surface.centroid))^2-norm(ray.origin-surface.centroid)^2-surface.extent^2>0 && return  nullIntersection,0 #Check if the ray goes anywhere near the sphere
    nIntersect = 0
    closestIntersection = nullIntersection
    for triangle in surface.Triangles
      intersection = polygonIntersect(ray,triangle)
      if intersection.isHit ≠ true
        continue
      end
      nIntersect += 1

      if (lIntersectMin == 0.0) | (intersection.distance < closestIntersection.distance)
        closestIntersection = intersection

      end
    end
  return closestIntersection, nIntersect
end




function load_ply_file(fileName::String)
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
  surface = Surface(nodeCoords, triIndices, nTriangles)



  totalSurfaceArea = sum([triangle.area for triangle in surface.Triangles])
  return surface
end

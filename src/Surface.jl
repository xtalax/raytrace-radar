using LinearAlgebra
include("raytracer-types.jl")
function SurfaceIntersect(surface::Surface, ray::Ray) # Outputs the index of the intersected triangle and the point of intersection.
    lIntersectMin = 0.0
    for triangle in surface.Triangles
      intersection = triangleIntersect(ray,triangle)
      if intersection.isHit ≠ true
        continue
      end
      nIntersect += 1
      if (lIntersectMin == 0.0) | (lIntersect < lIntersectMin)
        lIntersectMin = intersection.distance
        intersectedTriangle = triangle

      end
    end
  return triangle, intersection, nIntersect
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
      nodeCoords[2,i-iHeader] = parse(Float64,xyz[2])
      nodeCoords[3,i-iHeader] = parse(Float64,xyz[3])

    elseif i > iHeader+nNodes
      ijk=xyz=collect(m.match for m in eachmatch(r"(\d+)", line))

      triIndices[1,i-iHeader-nNodes] = parse(Int, ijk[2])+1
      triIndices[2,i-iHeader-nNodes] = parse(Int, ijk[3])+1
      triIndices[3,i-iHeader-nNodes] = parse(Int, ijk[4])+1
    end
    i += 1
  end
  close(iFile)
  println("~~~ Done ~~~ \n\n~~~ Building Surface ~~~")
  println(nodeCoords[:,triIndices[1,1]],typeof(nodeCoords[:,triIndices[1,2]]))

  surface = Surface(nodeCoords, triIndices, nTriangles)



  totalSurfaceArea = sum([triangle.area for triangle in surface.Triangles])
  return nTriangles, surface, totalSurfaceArea
end


################################################################################
# main start
################################################################################
# load the ply file and build triangles.
fileName = "/home/zander/Documents/Julia/raytrace-radar/cow.ply"
nTriangles, surface = load_ply_file(fileName)

#starting point of the line of sight
pStart = [0.0, 10.0, 0.0]

#direction of line of sight, must be vector of length 1.
r = [0.0, -1.0, 0.0]

# return index of triangle that contains the intersection point and
# the coordinates of the intersection point.
# if no intersection is found, return iTriangleIntersect = -1 and [0,0,0]
#iIntersectedTriangle, pIntersec, lIntersect, nIntersect= SurfaceIntersect(allTriangles, pStart, r)

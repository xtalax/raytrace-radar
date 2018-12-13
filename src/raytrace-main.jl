using LinearAlgebra
using StaticArrays
using Rotations
using FFTW
using JLD
#using Plots
using Distributed
@everywhere using SharedArrays
using ParallelDataTransfer
#using Constants
#using PhysicalConstants

include("raytracertypes.jl")
include("Surface.jl")
include("TriangleIntersect.jl")
include("trace-refract.jl")


################################################################################
# main start
################################################################################
# load the ply file and build triangles.

path = "/home/zander/github/Julia/raytrace-radar"
fileName = path*"/src/cube.ply"
objectName = match(r"\/\w+\.ply", fileName)
RefrIndex = √2.3
cow = load_ply_file(fileName, RefrIndex)

# initialize the bundle
println(cow.extent)
bundle = Bundle(Point(-3.0,0.0,0.0), Point(1.0,0.0,0.0)) # using default values - see raytracer-types for more

# initialize the detector

detector = Detector(Point(-3.1,0.0,0.0),Point(1.0,0.0,0.0)) # using default values - see raytracer-types for more

# Trace

nTheta = 361
Rplus = Matrix(RotZ(2π/nTheta))
Rminus = Matrix(RotZ(-2π/nTheta))
println(Rplus^(nTheta-1))

freqs = range(75*10^9,stop=90*10^9,length=1501)
mcFreqData = Complex.(zeros(1501,nTheta,2))
mcRangeData = Complex.(zeros(1501,nTheta,2))


for t in 0:nTheta-1
    bundleOrigin = (Rplus^t)*[-3.0,0.0,0.0]
    bundleDirection =(Rplus^t)*[1.0,0.0,0.0]

    bundle = Bundle(Point(bundleOrigin), Point(bundleDirection)) # using default values - see raytracer-types for more
    detector = Detector(Point(bundleOrigin - bundleDirection/1000),Point(bundleDirection)) # using default values - see raytracer-types for more

    @time Trace!!(detector,bundle,cow) #Fuckin trace the shit out of it
    detections = length(detector.pixels[1].times)
    println("\n ---------- Trace completed for $t","°","! There were $detections detections! ---------- \n")
    @time @sync Sf,Pf = freqsweepdistributed(detector,freqs)
    print("\nFreqency sweep completed!\n")
    @time  St,Pt = [ifft(Sf), ifft(Pf)]
    println("IFT completed!")

    mcFreqData[:, t+1, 1] = Sf ; mcFreqData[:, t+1, 2] = Pf
    mcRangeData[:, t+1, 1] = St ; mcRangeData[:, t+1, 2] = Pt

    #plot(detector.pixels[1].times, detector.pixels[1].ampS)
    #break
end


save(path*"/Data/"*objectName.match[2:end-4]*".jld", "mcFreqData", mcFreqData, "mcRangeData", mcRangeData)
#=
plotly(size = (1501,nTheta))
heatmap(20*log10.(abs.(mcRangeData[:,:,1]))) # plots the TE simulated data

plotly(size = (1501,nTheta))
heatmap(20*log10.(abs.(mcRangeData[:,:,2]))) # plot teh TM data
plotly(size = (1501,nTheta-50))
heatmap(20*log10.(abs.(mcRangeData[:,50:end,2]))) # plot teh TM data
=#

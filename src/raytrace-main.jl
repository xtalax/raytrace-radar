using LinearAlgebra
using StaticArrays
using Rotations
using FFTW
using JLD
#using Plots
using Distributed
@everywhere using DistributedArrays
using ParallelDataTransfer

#using Constants
#using PhysicalConstants

include("raytracertypes.jl")
include("utils.jl")
include("Surface.jl")
include("TriangleIntersect.jl")
include("trace-refract.jl")

#############################################################################
# main start
#############################################################################
# load the ply file and build triangles.
path = "/home/zander/github/Julia/raytrace-radar"
fileName = path*"/src/cube.ply"
objectName = match(r"\/\w+\.ply", fileName)
RefrIndex = √2.3
nodecoords, triindices, ntriangles = load_ply_file(fileName)
const cow =  Surface(nodecoords, triindices, ntriangles, RefrIndex)



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

TimesVec = [[0.0]]
AmpSVec = [[0.0]]
AmpPVec= [[0.0]]

###########################################################################
# Loop over θ
##########################################################################

#@time for t in 0:nTheta-350
    bundleOrigin = (Rplus^t)*[-3.0,0.0,0.0]
    bundleDirection =(Rplus^t)*[1.0,0.0,0.0]

    bundle = Bundle(Point(bundleOrigin), Point(bundleDirection)) # using default values - see raytracer-types for more
    detector = Detector(Point(bundleOrigin - bundleDirection/1000),Point(bundleDirection)) # using default values - see raytracer-types for more

    Trace!!(detector,bundle,cow) #Fuckin trace the shit out of it
    detections = length(detector.pixels[1].times)
    println("\n ---------- Trace completed for $t","°","! There were $detections detections! ---------- \n")


    @time Sf,Pf = freqsweepdistributed(detector,freqs)
    print("\nFreqency sweep completed!\n")
    @time  St,Pt = [ifft(Sf), ifft(Pf)]
    println("IFT completed!")

    mcFreqData[:, t+1, 1] = Sf ; mcFreqData[:, t+1, 2] = Pf
    mcRangeData[:, t+1, 1] = St ; mcRangeData[:, t+1, 2] = Pt

    timess, ampSs, ampPs = duplicateflatten(detector::Detector)


    if t == 0
        TimesVec[1] = timess
        AmpSVec[1] = ampSs
        AmpPVec[1] = ampPs
    else
        push!(TimesVec, timess)
        push!(AmpSVec, ampSs)
        push!(AmpPVec, ampPs)
    end
    #plot(detector.pixels[1].times, detector.pixels[1].ampS)
    #break

    if draw
        outerdiv() << initscene()


    end
#end




save(path*"/Data/"*objectName.match[2:end-4]*".jld", "mcFreqData", mcFreqData, "mcRangeData", mcRangeData)
#=

plotly(size = (1501,nTheta))
heatmap(20*log10.(abs.(mcRangeData[:,:,1]))) # plots the TE simulated data

plotly(size = (1501,nTheta))
heatmap(20*log10.(abs.(mcFreqData[:,:,1]))) # plots the TE simulated data

plotly(size = (1501,nTheta))
heatmap(20*log10.(abs.(mcRangeData[:,:,2]))) # plot teh TM data

plotly(size = (1501,nTheta-50))
heatmap(20*log10.(abs.(mcRangeData[:,50:end,2]))) # plot teh TM data
=#

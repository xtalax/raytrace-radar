function rₛ(θᵢ,θₜ,n₁,n₂)
    a = n₁*cos(θᵢ)
    b = n₂*cos(θₜ)
    out = (a-b)/(a+b)
    return out
     return abs(out)>1 ? out/norm(out) : out
 end


function rₚ(θᵢ,θₜ,n₁,n₂)
    a = n₁*cos(θₜ)
    b = n₂*cos(θᵢ)
    out = (a-b)/(a+b)
    return out
     return abs(out)>1 ? out/norm(out) : out
end




R₂(ϕ) = [cos(ϕ) -sin(ϕ) ; sin(ϕ) cos(ϕ)]

function Fresnel!(bundle::Bundle, ray::Ray, interaction::Intersection,n₁,n₂)
    pruned = true
    if ray.budget == 0
        print("💀")
        return
    end
    normal = interaction.normal
    d = ray.direction
    TE = ray.TEBasis
    TH = d×TE
    TEnew = normalize(normal×d)

    S, P = [ray.ampS, ray.ampP]
    ## Change polarisation basis and handle shit -- Splits the ray in to 4 to handle the different polarisations
    if abs(TE⋅TEnew) != 1 && !isnan(TEnew)
        THnew = d×TEnew
        ee = TE⋅TEnew*S
        eh = TH⋅TEnew*P
        he = TE⋅THnew*S
        hh = ee*P

        eray = Ray(ray.budget, ray.origin, ray.direction, TEnew, ray, ray.ξ, ee, he, ray.ϕs, ray.ϕs, ray.n )
        hray = Ray(ray.budget, ray.origin, ray.direction, TEnew, ray, ray.ξ, eh, hh,ray.ϕp, ray.ϕp, ray.n )

        Fresnel!(bundle, eray, interaction, n₁, n₂)
        Fresnel!(bundle, hray, interaction, n₁, n₂)
        return
    else
        TEnew = TE
    end


#    println(" ray direction ", ray.direction, " normal ", interaction.normal )
    θᵢ = vectorAngle(-ray.direction,normal)
    if θᵢ > π/2
        θᵢ= abs(θᵢ-π)

    elseif θᵢ < -π/2
        θᵢ= θᵢ+π
    end

    θₜ = asin(Complex(n₁*sin(θᵢ)/n₂))

    rs = rₛ(θᵢ,θₜ,n₁,n₂)
    rp = rₚ(θᵢ,θₜ,n₁,n₂)

    Rs = abs(rs)^2
    Rp = abs(rp)^2

    ϕs = normalize.(rs)
    ϕp = normalize.(rp)

    directionReflected = normalize(d-2*(d⋅normalize(normal))⋅normalize(normal))

    #Deal With reflection
    reflAmpS = ray.ampS*Rs
    reflAmpP = ray.ampP*Rp

    if (abs(reflAmpS)>0.0000001 )| (abs(reflAmpP)>0.0000001)

         Rs = (abs(Rs)> 1.0) ? 1.0 : Rs ; Rp = (abs(Rp)>1.0) ? 1.0 : Rp

        reflected = Ray(ray.budget-1, interaction.point+directionReflected*0.0001, directionReflected, TEnew, ray, ray.ξ+(interaction.distance+0.0001)*n₁ ,reflAmpS, reflAmpP, ray.ϕs*ϕs, ray.ϕp*ϕp, n₁ )
        pushfirst!(bundle.Rays, reflected) ; bundle.N += 1
        print("∠")
        isnan(reflected.ampS) | isnan(reflected.ampP) && throw("Nans have infected the rays")
        pruned = false

    else
    end

    if isreal(θₜ)
        Ts = 1-Rs
        Tp = 1-Rs
        if (abs(Ts)> 1.0) | (abs(Tp)>1.0)
            throw("Ts Greater than 1 at $Ts, Tp Greater than 1 at $Tp, Rs was $Rs, Rp was $Rp S was $S, P was $P, ϕs was $ϕp, ϕs was $ϕp")
        end
        Δ = 1-((n₁/n₂)^2)*(normal×d)⋅(normal×d)
        if Δ ≥ 0
            directionRefracted = normalize((n₁/n₂)*(normal×(-normal×d))-normal*√Δ)
            refrAmpS = ray.ampS*Ts
            refrAmpP = ray.ampP*Tp

        else
            refrAmpS = 0.0
            refrAmpP = 0.0
        end

        if (abs(refrAmpS)>0.0000001) | (abs(refrAmpP)>0.0000001)
            refracted = Ray(ray.budget-1, interaction.point+directionRefracted*0.0001, directionRefracted, TEnew, ray, ray.ξ + (interaction.distance)*n₁ + 0.0001*n₂,refrAmpS, refrAmpP,ray.ϕs,ray.ϕp, n₂ )
            pushfirst!(bundle.Rays, refracted) ; bundle.N += 1
            #println("refracted ray has [|ampP|, |ampS|] = ", [abs(refrAmpP), abs(refrAmpS)])
            isnan(refracted.ampS) | isnan(refracted.ampP) && throw("Nans have infected the rays")
            print("〽")
            pruned = false
        else
            print("×")

        end
    end


    pruned && print("×")



end

function Detect!(detector::Detector, ray::Ray)
    for pixel in detector.pixels
        detection = polygonIntersect(ray, pixel.shape)
        if detection.isHit
            d = normalize(ray.direction)
            normal = normalize(detection.normal)
            TEBasis = normalize(pixel.TEBasis)
            TMBasis = normalize(-normal×TEBasis)

            rayTE = normalize(ray.TEBasis)
            rayTM = normalize(d×ray.TEBasis)

            ampEE = ray.ampS*rayTE⋅TEBasis
            ampHH = ray.ampP*rayTM⋅TMBasis
            ampHE = ray.ampS*rayTE⋅TMBasis
            ampEH = ray.ampP*rayTM⋅TEBasis
            if abs(ampEE) > 1.0
                throw("amplitudes too big! EE $ampEE")
            end
        #    println("detection! --- ******** TE $TEBasis TM $TMBasis, rayTE $rayTE rayTM $rayTM, ampS $ampS ampP $ampP")
            print("*")
            push!(pixel.times, (ray.ξ+detection.distance*ray.n)/c₀)
            push!(pixel.times, (ray.ξ+detection.distance*ray.n)/c₀)

            push!(pixel.ampS, ampEE)
            push!(pixel.ampP, ampHH)
            push!(pixel.ampP, ampHE)
            push!(pixel.ampS, ampEH)
            push!(pixel.ϕs, ray.ϕs)
            push!(pixel.ϕp, ray.ϕp)
            push!(pixel.ϕp, ray.ϕs)
            push!(pixel.ϕs, ray.ϕp)
            break
        end
    end
end

function Trace!!(detector, bundle, surface)
    maxN = bundle.N
    while bundle.N > 0

        length(bundle.Rays) == 0 && break

        ray = pop!(bundle.Rays) ; bundle.N -=1

        interaction, nIntersect = SurfaceIntersect(cow,ray)

        # ~~~ DETECTION LOGIC GOES HERE ~~~ #
        detection = polygonIntersect(ray, detector.shape)
        if detection.isHit & (interaction.isHit & (detection.distance < interaction.distance) | !interaction.isHit)
             Detect!(detector, ray)
             continue
         end

        if interaction.isHit
            n₁, n₂ = nIntersect%2 == 0 ? (1, cow.RefrIndex) : (cow.RefrIndex, 1)
            Fresnel!(bundle, ray, interaction, n₁, n₂)
        else
            print("×")
        end
        maxN = maxN > bundle.N ? maxN : bundle.N
        #println(bundle.N, " rays left in bundle")
    end
    println("\n The largest number of rays was $maxN")

    print("\n\n")
    println()
end

@everywhere function getrange(n)
    tid = myid()-1
    nt = nprocs()-1
    d , r = divrem(n, nt)
    from = (tid - 1) * d + min(r, tid - 1) + 1
    to = from + d - 1 + (tid ≤ r ? 1 : 0)
    from:to
end

function freqsweepdistributed(detector::Detector,freqs::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}})
    sPol = SharedArray(Complex.(zeros(length(freqs),1)))
    pPol = SharedArray(Complex.(zeros(length(freqs),1)))
    pixel = detector.pixels[1]
    times = pixel.times
    ampS = pixel.ampS
    ampP = pixel.ampP
    ϕs = pixel.ϕs
    ϕp = pixel.ϕp
    sendto(workers(), freqs=freqs)
    sendto(workers(), times=times)
    sendto(workers(), ampS=ampS)
    sendto(workers(), ampP=ampP)
    sendto(workers(), ϕs=ϕs)
    sendto(workers(), ϕp=ϕp)

    @async @distributed for k in 1:nprocs()-1
         for j in getrange(1501)
            φ = exp.(im.*times.*2π.*freqs[j]) # phase delay
            sPol[j] = sum(ampS.*ϕs.*φ)
            pPol[j] = sum(ampP.*ϕp.*φ)
        end

    end
    return Array(sPol), Array(pPol)
end

function freqsweepserial(detector::Detector, freqs::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}})
    sPol = Complex.(zeros(length(freqs),1))
    pPol = Complex.(zeros(length(freqs),1))
    pixel = detector.pixels[1]
    times = pixel.times
    ampS = pixel.ampS
    ampP = pixel.ampP
    ϕs = pixel.ϕs
    ϕp = pixel.ϕp

    for j in 1:1501
        φ = exp.(im.*times.*2π.*freqs[j]) # phase delay
        sPol[j] = sum(ampS.*ϕs.*φ)
        pPol[j] = sum(ampP.*ϕp.*φ)
    end

    return sPol, pPol
end

function freqsweepserial(detector::Detector, freqs::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}})
    sPol = Complex.(zeros(length(freqs),1))
    pPol = Complex.(zeros(length(freqs),1))
    pixel = detector.pixels[1]
    times = pixel.times
    ampS = pixel.ampS
    ampP = pixel.ampP
    ϕs = pixel.ϕs
    ϕp = pixel.ϕp

    for j in 1:1501
        φ = exp.(im.*times.*2π.*freqs[j]) # phase delay
        sPol[j] = sum(ampS.*ϕs.*φ)
        pPol[j] = sum(ampP.*ϕp.*φ)
    end

    return sPol, pPol
end

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
    d = ray.direction
    normal = interaction.normal

    planeParallel =normalize(d×normal)
    #println("planeParallel",planeParallel)
    ϑ = vectorAngle(planeParallel, ray.TEBasis)

    #println("rotation angle", vectorAngle(planeParallel, ray.TEBasis))
    if !isnan(ϑ)
        R = R₂(ϑ)
        S, P =  R*[ray.ampS*ray.ϕs, ray.ampP*ray.ϕs]
        ray.ampS, ray.ampP, ray.ϕs, ray.ϕp = [abs(S), abs(P), normalize(S), normalize(P)]
        newTE = planeParallel
    else
        newTE = ray.TEBasis
    end
#    println(" ray direction ", ray.direction, " normal ", interaction.normal )
    θᵢ = vectorAngle(-ray.direction,interaction.normal)
    #θᵢ = θᵢ > π ? θᵢ-π : θᵢ
    θᵢ = abs(θᵢ) > π/2 ? θᵢ-π : θᵢ

    #println(θᵢ)

    θₜ = asin(Complex(n₂*sin(θᵢ)/n₁))
    rs = rₛ(θᵢ,θₜ,n₁,n₂)
    rp = rₚ(θᵢ,θₜ,n₁,n₂)
    Rs = abs(rs)^2
    Rp = abs(rp)^2
    ϕs = normalize.(rs)
    ϕp = normalize.(rp)
#println(Rp, " ", Rs)
    if isreal(θₜ)
        Ts = 1-abs(rs)^2
        Tp = 1-abs(rp)^2
        if (Ts> 1.0) | (Tp>1.0)
            throw("Ts Greater than 1 at $Ts, Tp Greater than 1 at $Tp, Rs was $Rs, Rp was $Rp S was $S, P was $P, ϕs was $ϕp, ϕs was $ϕp")
        end
        Δ = 1-(n₁/n₂)^2*(normal×d)⋅(normal×d)
        if Δ ≥ 0
            directionRefracted = (n₁/n₂)*(normal×(-normal×d))-normal*√Δ
            refrAmpS = ray.ampS*Ts
            refrAmpP = ray.ampP*Tp

        else
            refrAmpS = 0.0
            refrAmpP = 0.0
        end

        if (abs(refrAmpS)>0.01) | (abs(refrAmpP)>0.01)
            refracted = Ray(interaction.point+directionRefracted*0.001, directionRefracted, newTE, ray, ray.tof+interaction.distance/ray.c+0.001,refrAmpS, refrAmpP,ray.ϕs,ray.ϕp, c₀/n₂ )
            pushfirst!(bundle.Rays, refracted) ; bundle.N += 1
            #println("refracted ray has [|ampP|, |ampS|] = ", [abs(refrAmpP), abs(refrAmpS)])
            isnan(refracted.ampS) | isnan(refracted.ampP) && throw("Nans have infected the rays")
            print("〽")
            pruned = false
        else
            print("×")

        end
    end

    directionReflected = d-2*(d⋅normalize(normal))⋅normalize(normal)
    #println(directionReflected)

    #println(S," ",P)
    #println("Rs was $Rs, Rp was $Rp, θᵢ is $θᵢ, θₜ is $θₜ")

    #Deal With reflection
    reflAmpS = ray.ampS*Rs
    reflAmpP = ray.ampP*Rp

    if (abs(reflAmpS)>0.01 )| (abs(reflAmpP)>0.01)
        reflected = Ray(interaction.point+directionReflected*0.001, directionReflected, newTE, ray, ray.tof+interaction.distance/ray.c + 0.001,reflAmpS, reflAmpP, ray.ϕs*ϕs, ray.ϕp*ϕp, ray.c )
        pushfirst!(bundle.Rays, reflected) ; bundle.N += 1
        print("∠")
        isnan(reflected.ampS) | isnan(reflected.ampP) && throw("Nans have infected the rays")
        pruned = false

    else
    end
    pruned && print("×")



end

function Detect!(detector::Detector, ray::Ray)
    for pixel in detector.pixels
        detection = polygonIntersect(ray, pixel.shape)
        if detection.isHit
            d = ray.direction
            normal = detection.normal
            TEBasis = pixel.TEBasis
            TMBasis = normalize(normal×TEBasis)

            rayTE = ray.TEBasis
            rayTM =normalize(d×ray.TEBasis)

            ampS = (rayTE*ray.ampS*ray.ϕs + rayTM*ray.ampP*ray.ϕs)⋅TEBasis
            ampP = (rayTE*ray.ampS*ray.ϕs + rayTM*ray.ampP*ray.ϕp)⋅TMBasis
        #    println("detection! --- ******** TE $TEBasis TM $TMBasis, rayTE $rayTE rayTM $rayTM, ampS $ampS ampP $ampP")
            print("*")
            push!(pixel.times, ray.tof+detection.distance/ray.c)
            push!(pixel.ampS, abs(ampS))
            push!(pixel.ampP, abs(ampP))
            push!(pixel.ϕs, normalize(ampS) )
            push!(pixel.ϕp, normalize(ampP) )
            break
        end
    end
end

function Trace!!(detector, bundle, surface)
    while bundle.N > 0
        ray = pop!(bundle.Rays) ; bundle.N -=1
        # ~~~ DETECTION LOGIC GOES HERE ~~~ #
        detection = polygonIntersect(ray, detector.shape)
        if detection.isHit
             Detect!(detector, ray)
             #println("Detection! ray tof = ",ray.tof," Strengths = ",[ray.ampS,ray.ampP])
             continue
         end
        interaction, nIntersect = SurfaceIntersect(cow,ray)

        if interaction.isHit
            n₁, n₂ = nIntersect%2 == 1 ? (cow.RefrIndex, 1) :  (1, cow.RefrIndex)
            Fresnel!(bundle, ray, interaction, n₁, n₂)
        else
            print("×")
        end

        #println(bundle.N, " rays left in bundle")
    end
    print("\n\n")
end


function freqsweep(detector,freqs)
    sPol = Complex.(zeros(length(freqs),1))
    pPol = Complex.(zeros(length(freqs),1))
    for j in eachindex(freqs)
        pixel = detector.pixels[1]
        for i in eachindex(pixel.times)
            φ = exp(im*pixel.times[i]*2π*freqs[j]) # phase delay
            sPol[j] += pixel.ampS[i]*pixel.ϕs[i]*φ
            pPol[j] += pixel.ampP[i]*pixel.ϕp[i]*φ
        end
    end
    return sPol, pPol
end

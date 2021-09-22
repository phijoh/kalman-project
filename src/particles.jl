
function selectrandomparticles(frames, nparticles::Int64)

    T, width, height = size(frames)

    particles = zeros(Int64, nparticles, 2)
    particles[:,1] = sample(1:width, nparticles; replace=false)
    particles[:,2] = sample(1:height, nparticles; replace=false)

    weights = ones(nparticles)./nparticles

    return particles, weights

end


function getparticlevalues(frame::Matrix{Float64}, particles::Matrix{Int64}; rfsize=2)
    
    #TODO: smooth whole frame?

    width, height = size(frame)

    particlevalues = zeros(size(particles,1))

    for (i,particle) in eachrow(particles) |> enumerate

        x₀,y₀ = max.(particle .- (rfsize-1), 1)
        x₁,y₁ = min.(particle .+ (rfsize-1), width)

        particlevalues[i] = mean(frame[x₀:x₁,y₀:y₁])

    end

    return particlevalues

end


function duplicateparticles!(particles::Matrix{Int64}, weights::Vector{Float64})

    threshold = mean(weights)
    discardedparticleidx = weights .< threshold

    bestparticles = particles[.!discardedparticleidx,:]

    for i in (1:nparticles)[discardedparticleidx]
        
        duplicate = sample(1:size(bestparticles,1))
        particles[i,:] = bestparticles[duplicate,:]
    
    end

end


function updateweights(weights,probabilities)

end
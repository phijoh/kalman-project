
function selectrandomparticles(framesize, nparticles::Int64)

    width, height = framesize

    particles = zeros(Int64, nparticles, 4)
    x = sample(1:width, nparticles)
    y = sample(1:height, nparticles)

    u = sample(-(width ÷ 2):(width ÷ 2), nparticles)
    v = sample(-(height ÷ 2):(height ÷ 2), nparticles)

    particles[:, 1] = x
    particles[:, 2] = y
    particles[:, 3] = u
    particles[:, 4] = v

    weights = ones(nparticles) ./ nparticles

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
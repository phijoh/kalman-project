
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

"""
Get the luminance around (rfsize) the particles
"""
function getparticlevalues(frame::Matrix{Float64}, particles::Matrix{Int64}; rfsize=2)
    
    #TODO: smooth whole frame?

    width = size(frame, 1) # FIXME: What happens with height?

    particlevalues = zeros(size(particles,1))

    for (i, particle) in eachrow(particles) |> enumerate
        (x, y, u, v) = particle
        x₀,y₀ = max.([x, y] .- (rfsize-1), 1)
        x₁,y₁ = min.([x, y] .+ (rfsize-1), width)

        particlevalues[i] = mean(frame[x₀:x₁,y₀:y₁])

    end

    return particlevalues

end

"""
Drop particles with low weight and substitute with high weight particles
"""
function duplicateparticles!(particles::Matrix{Int64}, weights::Vector{Float64})

    N = size(particles, 1)

    threshold = mean(weights)
    discardedparticleidx = weights .< threshold

    bestparticles = particles[.!discardedparticleidx,:]

    for i in (1:N)[discardedparticleidx]
        
        duplicate = sample(1:size(bestparticles,1))
        particles[i,:] = bestparticles[duplicate,:]
    
    end

end

function moveparticles(particles::Matrix{Int64}, rest...)::Matrix{Int64}
    newparticles = copy(particles)
    moveparticles!(newparticles, rest...)

    return newparticles
end
function moveparticles!(particles::Matrix{Int64}, Σ::Matrix{Float64}, framesize::Tuple{Int64, Int64})

    width, height = framesize
    N = size(particles, 1)
    ν = rand(MvNormal(zeros(4), Σ), N)

    x = particles[:, 1:2]
    V = particles[:, 3:4]

    x′ = x + V + ν[1:2, :]'
    V′ = V + ν[3:4, :]'

    particles[:, 1] = @. round(Int64, mod(x′[:, 1], width))
    particles[:, 2] = @. round(Int64, mod(x′[:, 2], height))
    
    particles[:, 3:4] = @. round(Int64, V′)


end


"""
Update particles and weights given the frames, a current time t, a covariance matrix of motion Σ, and an observation noise.
"""
function step!(particles, w, frames, t, Σ, σ²ᵢ)
    frame = frames[t, :, :]
    frame′ = frames[t + 1, :, :]

    I = getparticlevalues(frame, particles)

    moveparticles!(particles, Σ, size(frame))

    I′ = getparticlevalues(frame′, particles)

    p = pdf.(Normal(0, σ²ᵢ), I - I′)

    w = normalize(w .* p)        

    duplicateparticles!(particles, w)
end

end
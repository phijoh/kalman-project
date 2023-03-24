"""
Estimates N particles using the data contained in the frames.
"""
function estimateparticles(
    N, frames::Vector{Frame}; 
    quantilecutoff = 0.9, rfsize = 5,
    verbose=false)

    T = length(frames)

    σ²ᵢ = mean(var.(frames)) # Get the luminance variance

    # Initialize weights and particles
    particles₀, weights₀ = samplerandomparticles(size(first(frames)), N)

    weightsovertime = Array{Float64}(undef, T, N)
    weightsovertime[1, :] = weights₀

    particlesovertime = Array{Int64}(undef, T, N, 4)
    particlesovertime[1, :, :] = particles₀

    # Initialize prior
    Σ = getvariance(particles₀, weights₀)

    for t in 2:T
        verbose && print("Iteration t = $t / $T \r")

        # One particle step with a random sample from the prior
        particles = particlesovertime[t-1, :, :]
        w = weightsovertime[t-1, :]
        fr, fr′ = frames[t-1], frames[t]

        particles, weights = particlestep(
            particles, w, fr, fr′;
            Σ, σ²ᵢ, rfsize, quantilecutoff
        )

        Σ .= getvariance(particles, weights)

        # Save particles and chain
        particlesovertime[t, :, :] = particles
        weightsovertime[t, :] = weights
    end

    return particlesovertime, weightsovertime
end


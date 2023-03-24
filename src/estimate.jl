"""
Estimates N particles using the data contained in the frames.
"""
function estimateparticles(
    N, frames::Vector{Frame}; 
    dimensions=4, quantilecutoff = 0.9, rfsize = 3,
    verbose=false)

    T = length(frames)

    σ²ᵢ = mean(var.(frames)) # Get the luminance variance

    # Initialize weights and particles
    particles₀, weights₀ = samplerandomparticles(size(first(frames)), N)

    weightsovertime = Array{Float64}(undef, T, N)
    weightsovertime[1, :] = weights₀

    particlesovertime = Array{Int64}(undef, T, N, dimensions)
    particlesovertime[1, :, :] = particles₀

    # Initialize prior
    Σ = getvariance(particles₀, weights₀)

    for tᵈ in 2:T
        verbose && print("Iteration tᵈ = $tᵈ / $T \r")

        # One particle step with a random sample from the prior
        particles = particlesovertime[tᵈ-1, :, :]
        w = weightsovertime[tᵈ-1, :]
        fr, fr′ = frames[tᵈ-1], frames[tᵈ]

        particles, weights = particlestep(
            particles, w, fr, fr′;
            Σ, σ²ᵢ, rfsize, quantilecutoff
        )

        Σ .= getvariance(particles, weights)

        # Save particles and chain
        particlesovertime[tᵈ, :, :] = particles
        weightsovertime[tᵈ, :] = weights
    end

    return particlesovertime, weightsovertime
end


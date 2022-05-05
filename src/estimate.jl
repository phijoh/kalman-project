function estimateparticle(T, N, frames; dimensions=4, verbose=false, rfsize)

    width, height = size(first(frames))
    σ²ᵢ = mean(var.(frames))

    # Initialize weights and particles
    particles₀, weights₀ = selectrandomparticles((width, height), N)

    weightsovertime = zeros(T, N)
    weightsovertime[1, :] = weights₀

    particlesovertime = zeros(Int64, T, N, dimensions)
    particlesovertime[1, :, :] = particles₀

    # Initialize prior
    Σ = diagm(ones(dimensions))

    for tᵈ in 2:T
        verbose && print("Iteration tᵈ = $tᵈ / $T \r")

        # One particle step with a random sample from the prior

        particles, weights = particlestep(
            particlesovertime[tᵈ-1, :, :], weightsovertime[tᵈ-1, :],
            frames[tᵈ-1], frames[tᵈ];
            Σ, σ²ᵢ, rfsize=rfsize
        )

        # TODO: Use maximum likelihood estimation

        Σ .= getvariance(particles, weights)

        # Save particles and chain
        particlesovertime[tᵈ, :, :] = particles
        weightsovertime[tᵈ, :] = weights
    end

    return particlesovertime, weightsovertime
end
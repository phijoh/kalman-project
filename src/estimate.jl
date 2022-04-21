function stepwiseparticle(T, N, frames; dimensions=4, verbose=false, rfsize)

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

    for t in 2:T
        verbose && print("Iteration t = $t / $T\r")

        # One particle step with a random sample from the prior

        particles, weights = particlestep(
            particlesovertime[t-1, :, :], weightsovertime[t-1, :],
            frames[t-1], frames[t];
            Σ, σ²ᵢ, rfsize=rfsize
        )

        # TODO: Use maximum likelihood estimation

        Σ .= getvariance(particles, weights)

        # Save particles and chain
        particlesovertime[t, :, :] = particles
        weightsovertime[t, :] = weights
    end

    return particlesovertime, weightsovertime
end
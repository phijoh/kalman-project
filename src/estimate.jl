@model function movement(p′, p, Σprior)

    N = size(p, 1)

    Σ ~ Σprior

    for i in 1:N
        p′[i, :] ~ MvNormal(p[i, :], Σ)
    end
    
end

function stepwiseparticle(T, N, frames; dimensions = 4, verbose = false, L = 100, sampler = NUTS(0.65))

    width, height = size(first(frames))

    chains = []

    # Initialize weights and particles
    particles₀, weights₀ = selectrandomparticles((width, height), N)
    
    weightsovertime = zeros(T, N)
    weightsovertime[1, :] = weights₀
    
    particlesovertime = zeros(Int64, T, N, dimensions)
    particlesovertime[1, :, :] = particles₀
    
    # Initialize prior
    Σ = MvLogNormal(MvNormal(zeros(dimensions), diagm(ones(dimensions))))

    for t in 2:T
        verbose && print("Iteration t = $t / $T")

        # One particle step with a random sample from the prior
        Σ̂ = rand(Σ) |> diagm

        particles, weights, losers = step(
            particlesovertime[t-1, :, :], weightsovertime[t-1, :], 
            frames, t, Σ̂, σ²ᵢ
        )

        # Estimate the posterior Σ using the winning particles
        notlosers = .~losers
        model = movement(particles[notlosers, :], particlesovertime[t-1, notlosers, :], Σ)

        chain = sample(model, sampler, L)

        # Make posterior new prior
        μₛ = zeros(dimensions)
        σ²ₛ = zeros(dimensions)

        for i in 1:dimensions
            μ = mean(chain.value[:, i])
            σ² = var(chain.value[:, i])

            μₛ[i], σ²ₛ[i] = logNtoN(μ, σ²)
        end

        Σ = MvLogNormal(MvNormal(μₛ , diagm(sqrt.(σ²ₛ))))

        # Save particles and chain
        particlesovertime[t, :, :] = particles
        weightsovertime[t, :] = weights
        push!(chains, chain)

    end

    return particlesovertime, weightsovertime, chains
end
function compensate(particles::Matrix{Int64}, τ::Int64, framesize::NTuple{2,Int64})

    compensatedparticles = copy(particles)
    dims = size(particles, 2)

    for _ in 1:τ
        moveparticles!(compensatedparticles, zeros(dims, dims), framesize)
    end

    return compensatedparticles

end

function sequencecompensation(
    particlesovertime::Array{Int64,3},
    weightsovertime::Matrix{Float64},
    τ::Int64, framesize::NTuple{2,Int64};
    verbose=false
)

    Tₑ, N, dims = size(particlesovertime)
    T = Tₑ + τ

    compensated = Array{Int64}(undef, T, N, dims)
    compweights = Matrix{Float64}(undef, T, N) # Weights associated with compensation

    # The first τ - 1 steps are simple movements of the first particle estimation
    compensated[1, :, :] .= particlesovertime[1, :, :]
    compweights[1, :] .= weightsovertime[1, :]

    for t ∈ 2:τ
        compensated[t, :, :] .= moveparticles(
            compensated[t-1, :, :],
            zeros(dims, dims), framesize
        )

        compweights[t, :] .= weightsovertime[1, :]
    end

    # After that there the tₑ particle is used to extrapolate to tₑ + τ
    for tₑ ∈ 1:Tₑ
        verbose && print("Compensating tₑ = $(tₑ + τ) / $(Tₑ + τ) \r")

        compensated[tₑ+τ, :, :] = compensate(
            particlesovertime[tₑ, :, :],
            τ, framesize
        )

        compweights[tₑ+τ, :, :] .= weightsovertime[tₑ, :, :]
    end

    return compensated, compweights

end
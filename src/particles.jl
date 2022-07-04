

function selectrandomparticles(framesize, nparticles::Int64)

    w, h = framesize

    particles = zeros(Int64, nparticles, 4)
    x = sample(1:w, nparticles)
    y = sample(1:h, nparticles)

    u = integernormal(nparticles, w ÷ 2)
    v = integernormal(nparticles, h ÷ 2)

    particles[:, 1] = x
    particles[:, 2] = y
    particles[:, 3] .= u
    particles[:, 4] .= v

    weights = ones(nparticles) ./ nparticles

    return particles, weights

end

"""
Get the luminance around (rfsize) the particles
"""
function getparticlevalues(frame::Frame, particles::Matrix{Int64}; rfsize=1)

    # TODO: smooth whole frame?

    width = size(frame, 1) # FIXME: What happens with height?
    N = size(particles, 1)

    particlevalues = Vector{Float64}(undef, N)

    for (i, particle) in eachrow(particles) |> enumerate
        x, y = particle[1:2]
        x₀, y₀ = max.([x, y] .- (rfsize - 1), 1)
        x₁, y₁ = min.([x, y] .+ (rfsize - 1), width)

        particlevalues[i] = mean(frame[x₀:x₁, y₀:y₁])

    end

    return particlevalues

end

function torus(x, w)
    mod(round(Int64, x) - 1, w) + 1
end

function moveparticles(particles::Matrix{Int64}, rest...)::Matrix{Int64}
    newparticles = copy(particles)
    moveparticles!(newparticles, rest...)

    return newparticles
end
function moveparticles!(particles::Matrix{Int64}, Σ::Matrix{Float64}, framesize::Tuple{Int64,Int64})

    width, height = framesize
    N = size(particles, 1)

    ν = all(Σ .≈ 0) ? zeros(4, N) : rand(MvNormal(zeros(4), sqrt.(Σ)), N)

    x = particles[:, 1:2]
    V = particles[:, 3:4]

    x′ = x + V + ν[1:2, :]'
    V′ = V + ν[3:4, :]'

    particles[:, 1] = torus.(x′[:, 1], width)
    particles[:, 2] = torus.(x′[:, 2], height)

    particles[:, 3:4] = @. round(Int64, V′)

end

Σ₀ = zeros(4, 4)

function likelihood(fr::Frame, fr′::Frame, particles, particles′, σ²ᵢ; rfsize)
    # Find likelihood of luminance
    I = getparticlevalues(fr, particles; rfsize)
    I′ = getparticlevalues(fr′, particles′; rfsize)
    return pdf.(Normal(0, σ²ᵢ), @. (I - I′))
end

"""
Update particles and weights given the current frame and the next frame, a covariance matrix of motion Σ, and the luminance noise.
"""
function particlestep(particles, w, fr::Frame, fr′::Frame; Σ, σ²ᵢ, intensity, rfsize)
    N = size(particles, 1)

    # Compute realized particles
    particles′ = moveparticles(particles, Σ, size(fr))
    L = likelihood(fr, fr′, particles, particles′, σ²ᵢ; rfsize)

    updated_weights = normalize(w .* L)
    cutoff = quantile(updated_weights, intensity)

    # Replace worst performing particles with theoretical particles
    losers = updated_weights .< cutoff
    winners = sample((1:N)[.!losers], sum(losers))

    particles′[losers, :] .= particles′[winners, :]
    updated_weights[losers] = updated_weights[winners]

    return particles′, normalize(updated_weights)
end

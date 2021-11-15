
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
function getparticlevalues(frame::Matrix{Float64}, particles::Matrix{Int64}; rfsize=10)
    
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

function torus(x, w)
    mod(round(Int64, x) - 1, w) + 1
end

function moveparticles(particles::Matrix{Int64}, rest...)::Matrix{Int64}
    newparticles = copy(particles)
    moveparticles!(newparticles, rest...)

    return newparticles
end
function moveparticles!(particles::Matrix{Int64}, Σ::Matrix{Float64}, framesize::Tuple{Int64, Int64})

    width, height = framesize
    N = size(particles, 1)

    ν = all(Σ .≈ 0) ? zeros(4, N) : rand(MvNormal(zeros(4), Σ), N)

    x = particles[:, 1:2]
    V = particles[:, 3:4]

    x′ = x + V + ν[1:2, :]'
    V′ = V + ν[3:4, :]'

    particles[:, 1] = torus.(x′[:, 1], width)
    particles[:, 2] = torus.(x′[:, 2], height)
    
    particles[:, 3:4] = @. round(Int64, V′)

end

Σ₀ = zeros(4, 4)

function likelihood(fr, fr′, particles, particles′, σ²ᵢ)
    # Find likelihood of luminance
    I = getparticlevalues(fr, particles)
    I′ = getparticlevalues(fr′, particles′)
    return pdf.(Normal(0, σ²ᵢ), @. (I - I′)^2)
end

"""
Update particles and weights given the frames, a current time t, a covariance matrix of motion Σ, and an observation noise.
"""
function step(particles, w, frames, t, Σ, σ²ᵢ)
    fr = frames[t, :, :]
    fr′ = frames[t + 1, :, :]

    # Compute theoretical particles
    particlesᴱ = moveparticles(particles, Σ₀, size(fr))

    p = likelihood(
        fr, fr′, 
        particles, particlesᴱ, σ²ᵢ)

    wᴱ = normalize(w .* p)

    # Compute realized particles
    particles′ = moveparticles(particles, Σ, size(fr))
    pₙ = likelihood(fr, fr′, particles, particles′, σ²ᵢ)

    wᴬ = normalize(w .* pₙ)        
    
    # Replace worst performing particles with theoretical particles
    losers = wᴬ .< wᴱ

    particles′[losers, :] .= particlesᴱ[losers, :]
    wᴬ[losers] = wᴱ[losers]

    return particles′, normalize(wᴬ)

end

filterframes = copy(frames)
for t in 1:T
    filterframes[t, :, :] = imfilter(frames[t, :, :], Kernel.LoG(25)) |> normalize
end

v₀ = 10.

Σ = [
    1 0 0 0;
    0 1 0 0;
    0 0 v₀ 0;
    0 0 0 v₀
]

T, width, height = size(frames)
particles, weights = selectrandomparticles((width, height), 2^16)
for t in 1:127
    println("Iteration t = $t / 127")
    global particles, weights = step(particles, weights, filterframes, t, Σ, σ²ᵢ)
end

plotparticles(particles, frames[127, :, :])

using JLD
using LinearAlgebra
using Random, Distributions, StatsBase
using Turing, ParticleFilters
using Parameters

using MAT, ImageFiltering
using Images, Rotations, CoordinateTransformations
using Base.Threads, Base.Iterators

using Printf
using Plots, StatsPlots
using DotEnv

# Environment
include("src/utilities/env.jl")
include("src/loadenv.jl")

# Utilities
include("src/utilities/matrix.jl") 
include("src/utilities/datautils.jl")
include("src/frame.jl")
# Data
include("src/loadposition.jl")

# Estimation
include("src/filtering.jl")
include("src/estimate.jl")
include("src/particles.jl")

# Diagnostic
include("src/diagnostic/chain.jl")
include("src/diagnostic/extrapolate.jl")
include("src/diagnostic/mcextrapolate.jl")

# Plots
include("src/plots/extrapolations.jl")
include("src/plots/particles.jl")

Random.seed!(seed)

Δt = 1 / framespersecond

function run(datapath, T, B, speed, duration, opacity; 
    sampler=NUTS(1_000, 0.65), dynamic=false, plotpath=plotpath, verbose=false)

    makeframes = loadgeneratedframe(datapath)
    frames = makeframes(speed, duration, opacity; dynamic=dynamic)

    verbose && println("Estimating position...")

    verbose && println("...done!")

    return x, x̂mc, x̂det, chain

end

duration, T, B = 128, 16, 1_000
speed = 1.2
opacity = 1.
dynamic = true

makeframes = loadgeneratedframe(datapath)
frames = makeframes(speed, duration, opacity; dynamic=dynamic)

wedgelength = 200
N = 2^12
T = length(frames)
width, height = size(frames[1])

Twedge = T - 80

σ²ᵢ = var(last(frames))
v₀ = tan(deg2rad(speed)) * wedgelength

Σ = [
    1. 0 0 0;
    0 1. 0 0;
    0 0 v₀ 0;
    0 0 0 v₀
]

particles₀, weights₀ = selectrandomparticles((width, height), N)

weightsovertime = zeros(Twedge, N); weightsovertime[1, :] = weights₀

particlesovertime = zeros(Int64, Twedge, N, 4); particlesovertime[1, :, :] = particles₀


for t in 2:Twedge
    print("Iteration t = $t / $Twedge \r")
    particles, weights = step(
        particlesovertime[t-1, :, :], 
        weightsovertime[t-1, :], 
        frames, t, Σ, σ²ᵢ
    )
    particlesovertime[t, :, :] = particles
    weightsovertime[t, :] = weights
end
println()
@gif for t in 1:Twedge
    print("Giffing $t / $Twedge \r")
    scatterparticles(particlesovertime[t, :, :], frames[t])
end

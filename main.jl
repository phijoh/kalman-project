using Base.Threads, Base.Iterators

using LinearAlgebra
using Random
using Distributions
using StatsBase

using Images
using CoordinateTransformations

using Plots

# Environment
using DotEnv
include("src/utilities/env.jl")
include("src/loadenv.jl")

# Models
include("src/illusions/twinklegoes.jl")

# Utilities
include("src/utilities/algorithms.jl")
include("src/utilities/matrix.jl")
include("src/utilities/datautils.jl")
include("src/utilities/particles.jl")

# Data
include("src/frame.jl")

# Estimation
include("src/particles.jl")
include("src/compensate.jl")
include("src/estimate.jl")

# Plots
include("src/plots/particles.jl")
include("src/plots/frame.jl")

Random.seed!(seed)

function runestimation(
    model::TwinkleGoesParameters; 
    N, τ = 60., verbose=false, 
    kwargs...)

    verbose && println("Generating frames...")

    frames = generateframes(model)

    verbose && println("...estimating position...")

    Tframes = length(frames)
    τframes = mstoframes(τ)

    particlesovertime, weightsovertime = estimateparticles(
        N, frames[1:(Tframes - τframes)];
        verbose, kwargs...
    )

    verbose && println("..compensating for neural delay τ...")
    compensated, compweights = trackwithdelay(
        particlesovertime, weightsovertime,
        τframes, size(last(frames));
        verbose
    )

    verbose && println("...done!")

    return compensated, compweights, frames

end

τ = 60.
speed = 2 # pixels / ms
opacity = 1
inducerduration = 500 # Estimation ms 
noiseduration = 180 # Overshoot ms

modelstatic = TwinkleGoesParameters(speed, opacity, false, inducerduration, noiseduration)
modeldynamic = TwinkleGoesParameters(speed, opacity, true, inducerduration, noiseduration)

N = 2^14
quantilecutoff = 0.9  # the quantile used for cutoff
rfsize = 10

results = runestimation.([modelstatic, modeldynamic]; τ, N, rfsize, quantilecutoff, verbose)

# Plot time series
models = [modelstatic, modeldynamic]
plotposition(results, models; labels = ["Static", "Dynamic"])
plotvelocity(results, models; labels = ["Static", "Dynamic"])

# Gif of localisation
particlesovertime, weightsovertime, frames = first(results)
anim = @animate for (t, frame) ∈ enumerate(frames)
    particles = results[1][t, :, :]
    weights = weightsovertime[t, :]

    plotparticledensity(particles, weights, frame)
end

gif(anim, joinpath(plotpath, "localisation.gif") )
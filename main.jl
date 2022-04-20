using JLD
using LinearAlgebra
using Random, Distributions, StatsBase
using Turing, ParticleFilters
using Parameters

using MAT, ImageFiltering
using Images, Rotations, CoordinateTransformations
using Base.Threads, Base.Iterators

using Printf, LaTeXStrings
using Plots, StatsPlots
using DotEnv

# Environment
include("src/utilities/env.jl")
include("src/loadenv.jl")

# Utilities
include("src/utilities/matrix.jl") 
include("src/utilities/datautils.jl")
include("src/utilities/distributions.jl")
include("src/utilities/particles.jl")

# Data
include("src/frame.jl")
include("src/loadposition.jl")

# Estimation
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

function run(datapath, T, speed, duration; opacity = 1., dynamic=true, plotpath=plotpath, N=2^12, verbose = false, kwargs...)

    verbose && println("Generating frames...")

    makeframes = loadgeneratedframe(datapath)
    frames = makeframes(speed, duration, opacity; dynamic=dynamic)

    verbose && println("...estimating position...")

    particlesovertime, weightsovertime, chains = stepwiseparticle(T, N, frames; verbose = verbose, kwargs...)

    verbose && println("...done!")

    return particlesovertime, weightsovertime, chains, frames

end

duration = 32 # Estimation frames 
overshoot = 16 # Overshoot frames
T = duration + overshoot # Total time
τ = 1 # Neural delay in frames (60ms) TODO: implement this.

speeds = [0.4] #, 0.6, 0.8, 1., 1.2]
S = length(speeds)
dimensions = 4
N = 2^12
dynamic = true
opacity = 0.5

results = Dict(
    :particles => Array{Float64}(undef, S, T, N, 4),
    :weights => Array{Float64}(undef, S, T, N),
    :speeds => speeds
)

for (s, speed) in enumerate(speeds)

    particlesovertime, weightsovertime, chains, frames = run(
        datapath, T, speed, duration; 
        opacity = opacity, L = 100, N = N,
        dynamic = dynamic, verbose = verbose
    )

    results[:particles][s, :, :, :] = particlesovertime
    results[:weights][s, :, :] = weightsovertime
end

scatterparticles(particlesovertime[T, :, :], frames[T])


fig = plotvariance(results, duration; dpi = 250)
savefig("plots/varplot_comparison.png")

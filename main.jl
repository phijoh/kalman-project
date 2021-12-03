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
include("src/utilities/distributions.jl")

# Data
include("src/frame.jl")
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

function run(datapath, T, B, speed, duration, opacity; dynamic=false, plotpath=plotpath, N=2^12, verbose = false, kwargs...)

    verbose && println("Generating frames...")

    makeframes = loadgeneratedframe(datapath)
    frames = makeframes(speed, duration, opacity; dynamic=dynamic)

    verbose && println("...estimating position...")

    particlesovertime, weightsovertime, chains = stepwiseparticle(T, N, frames; verbose = verbose)

    verbose && println("...done!")

    return particlesovertime, weightsovertime, chains

end

duration, T, B = 128, 16, 1_000
speed = 1.2
opacity = 1.
dynamic = true

particlesovertime, weightsovertime, chains = run(datapath, T, B, speed, duration, opacity; dynamic = true, verbose = true)
using JLD
using LinearAlgebra
using Random, Distributions
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
include("src/utilities/torus.jl")

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
include("src/plots/plotparticles.jl")

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

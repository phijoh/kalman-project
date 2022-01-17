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

function run(datapath, T, speed, duration, opacity; dynamic=false, plotpath=plotpath, N=2^12, verbose = false, kwargs...)

    verbose && println("Generating frames...")

    makeframes = loadgeneratedframe(datapath)
    frames = makeframes(speed, duration, opacity; dynamic=dynamic)

    verbose && println("...estimating position...")

    particlesovertime, weightsovertime, chains = stepwiseparticle(T, N, frames; verbose = verbose, kwargs...)

    verbose && println("...done!")

    return particlesovertime, weightsovertime, chains, frames

end

duration = 32
T = duration + 16
opacity = 1.
dynamic = true

for speed ∈ [0.9, 1., 1.1]

    particlesovertime, weightsovertime, chains, frames = run(
        datapath, T, speed, duration, opacity; 
        L = 100,
        dynamic = dynamic, verbose = true)

    title = "Speed = $speed"

    fig = plotvariance(particlesovertime, weightsovertime, duration; title = title, dpi = 250)

    intspeed = replace(speed |> string, '.' => '_')
    savefig("figures/varplot_$(intspeed).png")

end

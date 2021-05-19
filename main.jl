using JLD
using LinearAlgebra
using Random, Distributions
using Turing, ParticleFilters
using Parameters

using MAT, ImageFiltering
using Images, Rotations, CoordinateTransformations
using Base.Threads

using Printf
using Plots, StatsPlots
using DotEnv

# Environment
include("src/utilities/env.jl")
include("src/loadenv.jl")

# Utilities
include("src/utilities/matrix.jl") 
include("src/utilities/datautils.jl")

# Data
include("src/loadposition.jl")

# Estimation
include("src/filtering.jl")
include("src/estimate.jl")

# Diagnostic
include("src/diagnostic/chain.jl")
include("src/diagnostic/extrapolate.jl")
include("src/diagnostic/mcextrapolate.jl")

# Plots
include("src/plots/extrapolations.jl")

Random.seed!(seed)

Δt = 1 / framespersecond

function run(datapath, T, B, speed, duration, opacity; 
    sampler=NUTS(1_000, 0.65), dynamic=false, plotpath=plotpath, verbose=false)

    makeframes = loadgeneratedframe(datapath)
    frames = makeframes(speed, duration, opacity; dynamic=dynamic)

    verbose && println("Estimating position...")

    x = filtering(frames[1:(duration - 1), :, :])
    u = getvelocity(x, Δt)

    verbose && println("Estimating and simulating model...")

    model = movement(x, u, Δt)

    chain = sample(model, sampler, 1000, verbose=verbose)
    
    x̂mc = mcextrapolate(x, u, chain, Δt; T=T, B=B)
    x̂det, _ = extrapolate(x, u, chain, Δt; T=T)

    verbose && println("...done!")

    return x, x̂mc, x̂det, chain

end


if shallplot
    Plots.scalefontsizes(0.6)

    duration, T, B = 128, 16, 1_000
    to10string = p -> @sprintf("%.0f", p * 10)

    for speed in [0.6, 1.2], opacity in [0.3, 0.8]

        params = "$(to10string(speed))_$(to10string(opacity))"

        x, x̂mc, x̂det, chain = run(datapath, T, B, speed, duration, opacity; dynamic=true)

        describechain(chain; verbose=verbose, plotpath=plotpath, filename="$(params)_chain")
    
        plotfirstlikelihood(x, x̂det, chain; plotpath=plotpath, filename="$(params)_like")
    
        plotmontecarlo(x, x̂mc; plotpath=plotpath, filename="$(params)_mc")
    
    end

    Plots.resetfontsizes()
end
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

# FIXME: Move function and plot name based on parameter
function simulate(speed, duration, opacity; dynamic=false, plotpath=plotpath, verbose=verbose)

    plotpath = joinpath(plotpath, speed < 1 ? "slow" : "fast")

    makeframes = loadgeneratedframe(datapath)
    frames = makeframes(speed, duration, opacity; dynamic=dynamic)

    x = filtering(frames[1:(duration - 1), :, :])

    u = getvelocity(x, Δt)

    verbose && println("Estimating model...")
    model = movement(x, u, Δt)

    chain = sample(model, NUTS(1000, 0.65), 1000, verbose=verbose)

    Plots.scalefontsizes(0.6)

    x̂mc = mcextrapolate(x, u, chain, Δt; T=16, B=1_000)
    x̂det, ûdet = extrapolate(x, u, chain, Δt; T=16)

    describechain(chain; verbose=verbose, plotpath=plotpath)
    
    plotfirstlikelihood(
        x, x̂det, chain; plotpath=plotpath, 
        xs=-0.5:0.01:0., ys=0.5:0.01:1.
    )
    
    plotmontecarlo(x, x̂mc; plotpath=plotpath)
    
    Plots.resetfontsizes()
    

    verbose && println("...done!")

    return x̂mc

end


x̂slow = simulate(0.6, 128, 1.)
x̂fast = simulate(1.2, 128, 1.)

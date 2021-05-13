using JLD
using LinearAlgebra
using Random, Distributions
using Turing, ParticleFilters

using MAT, ImageFiltering
using Base.Threads

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

Plots.scalefontsizes(0.6)

Random.seed!(seed)

Δt = 1 / framespersecond
seconscreen = 128 / 200 # FIXME: Figure out from .mat file
accrate = 0.95

Nframeswithwedge = ceil(Int64, inv(Δt) * seconscreen)

x = loadposition(
    datapath; 
    cache=cache, verbose=verbose, 
    framelimit=Nframeswithwedge)
    
u = getvelocity(x, Δt)

verbose && println("Estimating model...")
model = movement(x, u, Δt)

chain = sample(model, NUTS(accrate), 2000, verbose=verbose)

shallplot && describechain(chain; verbose=verbose, plotpath=plotpath)

verbose && println("...done!")

x̂, û = extrapolate(x, u, chain, Δt; T=5)

shallplot && plotfirstlikelihood(x̂, û, chain; plotpath=plotpath)

Plots.resetfontsizes()
using JLD
using LinearAlgebra
using Random, Distributions
using Turing, ParticleFilters

using MAT, ImageFiltering
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

Plots.scalefontsizes(0.6)

Random.seed!(seed)

Δt = 1 / framespersecond
seconscreen = 128 / 200 # FIXME: Figure out from .mat file

Nframeswithwedge = ceil(Int64, inv(Δt) * seconscreen) 

x = loadposition(
    datapath; 
    cache=cache, verbose=verbose, framelimit=Nframeswithwedge)
    
u = getvelocity(x, Δt)

verbose && println("Estimating model...")
model = movement(x, u, Δt)

chain = sample(model, NUTS(1000, 0.65), 1000, verbose=verbose)

if shallplot

    x̂mc = mcextrapolate(x, u, chain, Δt; T=10, B=1000)
    x̂det, ûdet = extrapolate(x, u, chain, Δt; T=5)

    describechain(chain; verbose=verbose, plotpath=plotpath)
    
    plotfirstlikelihood(
        x, x̂det, chain; plotpath=plotpath, 
        xs=0.:0.01:0.5, ys=0.:0.01:0.5
    )
    
    plotmontecarlo(x, x̂mc; plotpath=plotpath)
    
end

Plots.resetfontsizes()
verbose && println("...done!")
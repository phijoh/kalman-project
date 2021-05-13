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

Random.seed!(seed)

Δt = 1 / framespersecond
seconscreen = 128 / 200 # FIXME: Figure out from .mat file

Nframeswithwedge = ceil(Int64, inv(Δt) * seconscreen)

x = loadposition(
    datapath; 
    cache=cache, verbose=verbose, 
    framelimit=Nframeswithwedge)
    
u = getvelocity(x, Δt)

if shallplot
    verbose && println("Estimating model...")
    model = movement(x, u, Δt)

    chain = sample(model, NUTS(0.65), 2000, verbose=verbose)
        
    Dᵥest = mean(chain[:Dᵥ])
    Dₓest = mean(chain[:Dₓ])
    σₚ²est = mean(chain[:σₚ²])

    σᵥest = inv(inv(σₚ²est) + inv(Dᵥest))
    γest = inv(1 + Dᵥest / σₚ²est)

    println("σᵥest = $σᵥest")
    println("γest = $γest")

    if !isnothing(plotpath)

        plot(chain, dpi=200)
        
        filename = joinpath(plotpath, "estimation")
        savefig(filename)

    else println("No plot path provided in .env") end

    verbose && println("...done!")

end
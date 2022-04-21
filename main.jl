using JLD
using LinearAlgebra
using Random, Distributions, StatsBase

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

# Plots
include("src/plots/extrapolations.jl")
include("src/plots/particles.jl")

Random.seed!(seed)

function runestimation(datapath, T, speed, duration; opacity=1.0, dynamic=true, N=2^12, verbose=false, rfsize=1)

    verbose && println("Generating frames...")

    makeframes = loadgeneratedframe(datapath)
    frames = makeframes(speed, duration, opacity; dynamic=dynamic)

    verbose && println("...estimating position...")

    particlesovertime, weightsovertime = stepwiseparticle(T, N, frames; verbose=verbose, rfsize=rfsize)

    verbose && println("...done!")

    return particlesovertime, weightsovertime, frames

end

duration = 48 # Estimation frames 
overshoot = 32 # Overshoot frames
rfsize = 1
T = duration + overshoot # Total time
τ = 1 # Neural delay in frames (60ms) TODO: implement this.

specifications = [
    (1.2, 1.0, false),
    (1.2, 1.0, true)
] # speed, opacity, dynamic noise present

S = length(specifications)
dimensions = 4
N = 2^12

results = Dict(
    :particles => Array{Int64}(undef, S, T, N, dimensions),
    :weights => Array{Float64}(undef, S, T, N),
    :frames => Array{Frame}[],
    :specs => specifications,
)

for (s, specs) in enumerate(specifications)

    speed, opacity, dynamic = specs

    particlesovertime, weightsovertime, frames = runestimation(
        datapath, T, speed, duration;
        opacity=opacity, N=N,
        dynamic=dynamic, verbose=verbose,
        rfsize=rfsize
    )

    push!(results[:frames], frames)

    results[:particles][s, :, :, :] = particlesovertime
    results[:weights][s, :, :] = weightsovertime
end

if shallplot

    precfig = plotprecision(results, duration; dpi=250)
    savefig(precfig, joinpath(plotpath, "precision.png"))

    angleplot = plotangle(results, duration; after=5, dpi=250, labels=["static", "dynamic"], marker=:o, legend=:topleft)
    savefig(angleplot, joinpath(plotpath, "extrapolation.png"))



    for (s, specs) ∈ enumerate(specifications)
        speed, opacity, dynamic = specs
        specname = replace(join(specs, "-"), "." => "_")

        # Gif of first specification
        frames = results[:frames][s]
        particlesovertime = results[:particles][s, :, :, :]
        weightsovertime = results[:weights][s, :, :]

        anim = @animate for t ∈ 1:T
            plotexpectedposition(particlesovertime[t, :, :], weightsovertime[t, :], frames[t])
        end

        gif(anim, joinpath(plotpath, "estpos-$specname.gif"), fps=15)
    end

end
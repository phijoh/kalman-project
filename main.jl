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
include("src/compensate.jl")

# Plots
include("src/plots/extrapolations.jl")
include("src/plots/particles.jl")

Random.seed!(seed)

function runestimation(inducerduration, noiseduration; τ, speed, opacity, dynamic, N, verbose=false, rfsize=1, dimensions=4)

    verbose && println("Generating frames...")

    frames = makeframes(inducerduration, noiseduration, speed, opacity; dynamic=dynamic)

    verbose && println("...estimating position...")

    T = length(frames)
    particlesovertime, weightsovertime = estimateparticle(
        T - τ, N, frames;
        dimensions, verbose=verbose, rfsize=rfsize
    )

    verbose && println("..compensating for neural delay τ...")
    compensated, compweights = sequencecompensation(
        particlesovertime, weightsovertime,
        τ, size(last(frames))
    )

    verbose && println("...done!")

    return compensated, compweights, frames

end

inducerduration = 640 # Estimation ms 
noiseduration = 100 # Overshoot ms

Tᵢ = mstoframes(inducerduration)
Tₙ = mstoframes(noiseduration)
makeframes = loadgeneratedframe(datapath)

T = Tᵢ + Tₙ # Total time
τₘ = mstoframes(100) # Neural delay in ms, TODO: implement this.
rfsize = 1

specifications = [
    (0.16, 1.0, false, τₘ)
# (0.16, 1.0, true)
] # speed, opacity, dynamic noise, τ

S = length(specifications)
dimensions = 4
N = 2^15

results = Dict(
    :particles => Array{Int64}(undef, S, T, N, dimensions),
    :weights => Array{Float64}(undef, S, T, N),
    :frames => Array{Frame}[],
    :specs => specifications,
)

for (s, specs) in enumerate(specifications)

    speed, opacity, dynamic, τ = specs

    compensatedparticles, weightsovertime, frames = runestimation(
        inducerduration, noiseduration;
        τ, speed, opacity,
        N, dynamic, verbose,
        rfsize
    )

    push!(results[:frames], frames)

    results[:particles][s, :, :, :] = compensatedparticles
    results[:weights][s, :, :] = weightsovertime
end



if shallplot
    precfig = plotprecision(results, Tᵢ)
    savefig(precfig, joinpath(plotpath, "precision.png"))

    angleplot = plotangle(results, Tᵢ; after=5, labels=["no delay", "short delay", "delay"], marker=:o, legend=:topleft)
    savefig(angleplot, joinpath(plotpath, "extrapolation.png"))


    for (s, specs) ∈ enumerate(specifications)
        speed, opacity, dynamic, τ = specs
        specname = replace(join(specs, "-"), "." => "_")

        # Gif of first specification
        frames = results[:frames][s]
        particlesovertime = results[:particles][s, :, :, :]
        weightsovertime = results[:weights][s, :, :]

        anim = @animate for t ∈ 1:T
            particles = particlesovertime[t, :, :]
            weights = weightsovertime[t, :]

            plotexpectedposition(particles, weights, frames[t])
        end

        gif(anim, joinpath(plotpath, "estpos-$specname.gif"), fps=15)
    end

end

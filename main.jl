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
include("src/utilities/algos.jl")
include("src/utilities/matrix.jl")
include("src/utilities/datautils.jl")
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

function runestimation(inducerduration, noiseduration; τ, speed, opacity, dynamic, N, verbose=false, rfsize=1, dimensions=4, intensity=0.5)

    verbose && println("Generating frames...")

    frames = makeframes(inducerduration, noiseduration, speed, opacity; dynamic=dynamic)

    verbose && println("...estimating position...")

    T = length(frames)
    particlesovertime, weightsovertime = estimateparticle(
        T - τ, N, frames;
        dimensions, verbose, rfsize, intensity
    )

    verbose && println("..compensating for neural delay τ...")
    compensated, compweights = sequencecompensation(
        particlesovertime, weightsovertime,
        τ, size(last(frames));
        verbose
    )

    verbose && println("...done!")

    return compensated, compweights, frames

end

inducerduration = 900 # Estimation ms 
noiseduration = 400 # Overshoot ms

Tᵢ = ceil(Int64, mstoframes(inducerduration))
Tₙ = ceil(Int64, mstoframes(noiseduration))
makeframes = loadgeneratedframe(datapath)

T = Tᵢ + Tₙ # Total time
τₘ(ms) = ceil(Int64, mstoframes(ms)) # Neural delay in ms, TODO: implement this.
rfsize = 1

specifications = product(
                     [0.16], # Speeds
                     [1.0], # Opacity
                     [false, true], # Dynamic noise
                     [τₘ(60)] # Neural delay
                 ) |> collect |> vec # Necessary to preserve order

S = length(specifications)
dimensions = 4
N = 2^16
intensity = 0.95  # the quantile used for cutoff

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
        rfsize, intensity
    )

    push!(results[:frames], frames)

    results[:particles][s, :, :, :] = compensatedparticles
    results[:weights][s, :, :] = weightsovertime
end

verbose && println("Plotting...")
precfig = plotprecision(results, Tᵢ; dpi=180)
savefig(precfig, joinpath(plotpath, "precision.png"))

anglefig = plotangle(results, Tᵢ; dpi=180)
savefig(anglefig, joinpath(plotpath, "angle.png"))


if shallplot


    for (s, specs) ∈ enumerate(specifications)

        verbose && println("Making gif for specification $(s) / $(length(specifications))...")

        speed, opacity, dynamic, τ = specs
        specname = replace(join(specs, "-"), "." => "_")

        frames = results[:frames][s]
        particlesovertime = results[:particles][s, :, :, :]
        weightsovertime = results[:weights][s, :, :]

        anim = @animate for t ∈ 1:T
            particles = particlesovertime[t, :, :]
            weights = weightsovertime[t, :]

            plotexpectedposition(particles, weights, frames[t]; title=latexstring("\$ t = $(t - τ) \$"))
        end

        gif(anim, joinpath(plotpath, "estpos-$specname.gif"), fps=15)
    end

end

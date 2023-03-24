using Base.Threads, Base.Iterators

using LinearAlgebra
using Random
using Distributions
using StatsBase

using Images
using CoordinateTransformations

using Plots

# Environment
using DotEnv
include("src/utilities/env.jl")
include("src/loadenv.jl")

# Models
include("src/illusions/twinklegoes.jl")

# Utilities
include("src/utilities/algorithms.jl")
include("src/utilities/matrix.jl")
include("src/utilities/datautils.jl")
include("src/utilities/particles.jl")

# Data
include("src/frame.jl")

# Estimation
include("src/particles.jl")
include("src/compensate.jl")
include("src/estimate.jl")

# Plots
include("src/plots/particles.jl")
include("src/plots/frame.jl")

Random.seed!(seed)

function runestimation(model::TwinkleGoesParameters; N, verbose=false, rfsize=1, dimensions=4, intensity=0.5)

    verbose && println("Generating frames...")

    frames = generateframes(model)

    verbose && println("...estimating position...")

    T = length(frames)
    particlesovertime, weightsovertime = estimateparticles(
        N, frames[1:(T - τ)];
        dimensions, verbose, rfsize, intensity
    )

    verbose && println("..compensating for neural delay τ...")
    compensated, compweights = trackwithdelay(
        particlesovertime, weightsovertime,
        τ, size(last(frames));
        verbose
    )

    verbose && println("...done!")

    return compensated, compweights, frames

end

inducerduration = 500 # Estimation ms 
noiseduration = 300 # Overshoot ms
model = TwinkleGoesParameters(2., 1., true, 500, 300)

Tᵢ = mstoframes(inducerduration)
Tₙ = mstoframes(noiseduration)

T = Tᵢ + Tₙ # Total time
rfsize = 5


specifications = product(
                     [2], # Speeds
                     [1.0], # Opacity
                     [false, true], # Dynamic noise
                     [mstoframes(60)] # Neural delay
                 ) |> collect |> vec # Necessary to preserve order

S = length(specifications)
dimensions = 4
N = 2^14
intensity = 0.9  # the quantile used for cutoff

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



# Plotting

verbose && println("Plotting...")
precfig = plotprecision(results, Tᵢ; dpi=180, legend=:topleft)
savefig(precfig, joinpath(plotpath, "precision.png"))

posfig = plotposition(results, Tᵢ, 2; dpi=180, legend=:topright, ylabel="\$x\$")
savefig(posfig, joinpath(plotpath, "position.png"))

velfig = plotposition(results, Tᵢ, 4; dpi=180, legend=:topleft, ylabel="\$v\$")
savefig(velfig, joinpath(plotpath, "velocity.png"))

# Make gif for debugging purposes
for (s, specs) ∈ enumerate(specifications)

    verbose && println("Making gif for specification $(s) / $(length(specifications))...")

    speed, opacity, dynamic, τ = specs
    specname = replace(join(specs, "-"), "." => "_")

    frames = results[:frames][s]
    particlesovertime = results[:particles][s, :, :, :]
    weightsovertime = results[:weights][s, :, :]

    anim = @animate for t ∈ 1:T
        p = particlesovertime[t, :, :]
        w = weightsovertime[t, :]
        fr = frames[t]

        plotparticledensity(p, w, fr; title="\$t = $(t - τ) \$")
    end

    gif(anim, joinpath(plotpath, "estpos-$specname.gif"), fps=15)
end

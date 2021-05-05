using Random, Distributions
using Turing
using StatsPlots

Random.seed!(1123)

include("src/dgp.jl")
include("src/estimate.jl")

σₓ, σᵧ = 8., 1.
vel₀ = zeros(2)

vel = dgp(σₓ, σᵧ)
model = movement(nothing, vel, vel₀, 0.05)
chain = sample(model, HMC(0.05, 10), 1000, progress=true)

p_summary = chain[:h]
plot(p_summary, seriestype=:histogram)
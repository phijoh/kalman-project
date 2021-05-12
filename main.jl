using LinearAlgebra
using Random, Distributions
using Turing, ParticleFilters

using Plots, StatsPlots

Random.seed!(1123)

include("src/utils.jl")
include("src/estimate.jl")
include("src/dgp.jl")

Dₓ = 0.2
Dᵥ = 0.25
σₚ² = 0.25

Δt = 0.05
x₀ = [100., 0.]
u₀ = [0., 0.]

x, u = dgp(x₀, u₀, Dₓ, Dᵥ, σₚ², Δt; T=4_000)

model = movement(x, u, Δt)
chain = sample(model, NUTS(0.65), 3_000)

Dᵥest = mean(chain[:Dᵥ])
Dₓest = mean(chain[:Dₓ])
σₚ²est = mean(chain[:σₚ²])

σᵥ = inv(inv(σₚ²) + inv(Dᵥ))
γ = inv(1 + Dᵥ / σₚ²)

σᵥest = inv(inv(σₚ²est) + inv(Dᵥest))
γest = inv(1 + Dᵥest / σₚ²est)

println("σᵥ = $σᵥ, σᵥest = $σᵥest")
println("γ = $γ, γest = $γest")


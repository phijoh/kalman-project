include("dgp.jl")
include("kalmanfilter.jl")

using Distributions, Random
using LinearAlgebra

using Plots

Random.seed!(121212)

pixels, T = 400, 60

σ = 5
Σ = [
    σ^2 10;
    10 σ^2
] # Matrix(σ^2 * I, 2, 2) # Variance covariance matrix

ε = MultivariateNormal(zeros(2), Σ)

Y = generatedata(pixels, T, ε)

Ŷ, P, posterror = kalmanfilter(Y; a₁=[0, 100], Σ=Σ)

plot(Y, c=:blue, label="true")
plot!(Ŷ, c=:red, label="estimate")
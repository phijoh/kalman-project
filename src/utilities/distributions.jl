function logNtoN(M, V)
    σ² = log(1 + V / (M^2))
    μ = log(M) - (σ² / 2)

    return μ, σ²
end
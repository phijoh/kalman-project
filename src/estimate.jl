@model function movement(x, u, Δt)

    T, M = size(x)

    # Parameters to estimate
    Dₓ ~ InverseGamma(2, 0.5)
    Dᵥ ~ InverseGamma(2, 0.5)
    σₚ² ~ InverseGamma(2, 0.5)

    for t = 1:T
        
        νₓ = Dₓ 
        νᵤ = inv(inv(σₚ²) + inv(Dᵥ))
        γ = inv(1 + Dᵥ / σₚ²)

        xₚ = t > 1 ? x[t - 1, :] : x[1, :]
        uₚ = t > 1 ? u[t - 1, :] : u[1, :]

        x[t, :] ~ MvNormal(xₚ + Δt * uₚ, νₓ * Δt)
        u[t, :] ~ MvNormal(γ * uₚ, νᵤ * Δt)

    end
    
end

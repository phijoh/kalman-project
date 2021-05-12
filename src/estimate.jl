@model movement(x, u, Δt) = begin

    T, M = size(x)

    # Parameters to estimate
    Dₓ ~ InverseGamma(2, 0.5)
    Dᵥ ~ InverseGamma(2, 0.5)
    σₚ² ~ InverseGamma(2, 0.5)

    νₓ = Dₓ 
    νᵤ = inv(inv(σₚ²) + inv(Dᵥ))
    γ = inv(1 + Dᵥ / σₚ²)

    x[1, :] ~ MvNormal(x[1, :], νₓ * Δt)
    u[1, :] ~ MvNormal(u[1, :], νᵤ * Δt)

    for t = 2:T
        
        νₓ = Dₓ 
        νᵤ = inv(inv(σₚ²) + inv(Dᵥ))
        γ = inv(1 + Dᵥ / σₚ²)

        xₚ = x[t - 1, :]
        uₚ = u[t - 1, :]

        x[t, :] ~ MvNormal(xₚ + Δt * uₚ, νₓ * Δt)
        u[t, :] ~ MvNormal(γ * uₚ , νᵤ * Δt)

    end
    
end


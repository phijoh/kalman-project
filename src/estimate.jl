@model function movement(x, u, Δt)

    T = size(x, 1)

    # Parameters to estimate
    Dₓ ~ InverseGamma(2, 1)
    Dᵥ ~ InverseGamma(2, 1)
    σₚ² ~ InverseGamma(2, 1)

    # Initial position
    x[1, :] ~ MvNormal([0., 0.], 1.)
    u[1, :] ~ MvNormal([0., 0.], 1.)

    for t = 2:T

        νₓ = Dₓ 
        νᵤ = inv(inv(σₚ²) + inv(Dᵥ))
        γ = inv(1 + Dᵥ / σₚ²)

        uₑ = u[t - 1, :] * γ
        xₑ = x[t - 1, :] + Δt * uₑ

        x[t, :] ~ MvNormal(xₑ, νₓ * Δt)
        u[t, :] ~ MvNormal(uₑ, νᵤ * Δt)

    end
    
end

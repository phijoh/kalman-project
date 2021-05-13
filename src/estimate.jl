@model function movement(x, u, Δt)

    T = size(x, 1)

    # Parameters to estimate
    Dₓ ~ InverseGamma(3, 300)
    Dᵥ ~ InverseGamma(3, 300)
    σₚ² ~ InverseGamma(3, 300)

    # Initial position
    x[1, :] ~ MvNormal([0., 0.], 2.)
    u[1, :] ~ MvNormal([0., 0.], 2.)

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

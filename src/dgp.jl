"""
Generates data with T time points
"""
function dgp(x₀, u₀, Dₓ, Dᵥ, σₚ², Δt; T=200)
    
    σᵥ = inv(σₚ²) + inv(Dᵥ)
    γ = inv(1 + Dᵥ / σₚ²)

    νₓ = MvNormal(zeros(2), Dₓ * Δt)
    νᵤ = MvNormal(zeros(2), inv(σᵥ) * Δt)

    x = zeros(T, 2); x[1, :] = x₀
    u = zeros(T, 2); u[1, :] = u₀

    for t in 2:T
        x[t, :] = x[t - 1, :] + u[t - 1, :] + rand(νₓ)
        u[t, :] = γ * u[t - 1, :] + rand(νᵤ)
    end

    return x, u
end


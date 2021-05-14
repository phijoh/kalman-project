function constructextrapolation(chain::Chains, Δt)
    Dᵥ = mean(chain, :Dᵥ)
    Dₓ = mean(chain, :Dₓ)
    σₚ² = mean(chain, :σₚ²)
    
    γ = inv(1 + Dᵥ / σₚ²)

    function step(x, u)

        x′ = x + Δt * u
        u′ = γ * u

        return x′, u′
    end
end

function extrapolate(x, u, chain, Δt; T=15)

    g = constructextrapolation(chain, Δt)

    x₀ = x[end, :]
    u₀ = u[end - 1, :]

    xforecast = zeros(T + 1, 2); xforecast[1, :] = x₀
    uforecast = zeros(T + 1, 2); uforecast[1, :] = u₀
    
    for t in 1:T

        xc = xforecast[t, :]
        uc = uforecast[t, :]
        
        x′, u′ = g(xc, uc)
        
        xforecast[t + 1, :] = x′
        uforecast[t + 1, :] = u′
        
    end

    x̂ = xforecast[2:end, :]
    û = uforecast[2:end, :]

    return x̂, û
end
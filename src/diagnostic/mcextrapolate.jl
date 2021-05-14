function constructmcextrapolation(Dᵥ, Dₓ, σₚ², Δt)
    γ = inv(1 + Dᵥ / σₚ²)
    νᵤ = inv(inv(σₚ²) + inv(Dᵥ))


    function step(x, u)

        x′ = rand(MvNormal(x + Δt * u, Dₓ * Δt))
        u′ = rand(MvNormal(γ * u, νᵤ * Δt))

        return x′, u′
    end 
end

function mcextrapolate(x, u, chain, Δt; T=15, B=100)
    
    mcx = zeros(2, T, B)
    mcu = copy(mcx)

    @threads for b in 1:B

        σₚ² = sample(chain[:σₚ²])
        Dₓ = sample(chain[:Dₓ])
        Dᵥ = sample(chain[:Dᵥ])

        g = constructmcextrapolation(Dᵥ, Dₓ, σₚ², Δt)

        for t in 1:T

            xprev = t > 1 ? mcx[:, t - 1, b] : x[end, :]
            uprev = t > 1 ? mcu[:, t - 1, b] : u[end - 1, :]
            
            x′, u′ = g(xprev, uprev)
            
            mcx[:, t, b] = x′
            mcu[:, t, b] = u′
            
        end
    end

    return mcx

end


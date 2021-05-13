function constructextrapolation(chain, Δt)
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

function plotfirstlikelihood(x̂, û, chain; plotpath="")
    
    likelihood = (x₁, x₂) -> pdf(MvNormal(x̂[1, :], mean(chain, :Dₓ)), [x₁, x₂])

    xs = -0.5:0.005:-0.1
    ys = reverse(-xs)
    
    heatmap(
        xs, ys, likelihood, 
        title="Position extrapolation with first point likelihood",
        xlabel="x", ylabel="y", legend=:none, aspect_ratio=1,
        xlims=extrema(xs), ylims=extrema(ys), dpi=300)

    scatter!(
        x[:, 1], x[:, 2], 
        label="estimate", legend=:bottomright,
        c=:white, markersize=2)

    scatter!(
        [x̂[1, 1]], [x̂[1, 2]],
        label="extrapolation", c=:black,
        markersize=2)

    savefig(joinpath(plotpath, "extrapolation"))

end
function describechain(c::Chains; verbose=false, plotpath=nothing)
    Dᵥ = mean(c, :Dᵥ)
    Dₓ = mean(c, :Dₓ)
    σₚ² = mean(c, :σₚ²)

    σᵥ = inv(inv(σₚ²) + inv(Dᵥ))
    γ = inv(1 + Dᵥ / σₚ²)

    if verbose

        println("
            Estimated Dᵥ = $(@sprintf("%.2f", Dᵥ)), Dₓ = $(@sprintf("%.2f", Dₓ)), σₚ² = $(@sprintf("%.2f", σₚ²))
            Derived σᵥ = $(@sprintf("%.2f", σᵥ)), γ = $(@sprintf("%.2f", γ))
        ")
    
    end

    if !isnothing(plotpath)

        plot(chain, dpi=200)
        
        filename = joinpath(plotpath, "estimation")
        savefig(filename)
    end
end

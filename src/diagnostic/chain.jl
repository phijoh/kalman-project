function describechain(c::Chains; verbose=false, plotpath=nothing)
    Dᵥ = mean(c, :Dᵥ)
    Dₓ = mean(c, :Dₓ)
    σₚ² = mean(c, :σₚ²)

    σᵥ = inv(inv(σₚ²) + inv(Dᵥ))
    γ = inv(1 + Dᵥ / σₚ²)

    verbose && println("Estimated σᵥ = $σᵥ, γ = $γ")

    if !isnothing(plotpath)

        plot(chain, dpi=200)
        
        filename = joinpath(plotpath, "estimation")
        savefig(filename)
    end
end
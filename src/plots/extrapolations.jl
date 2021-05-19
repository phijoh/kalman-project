function plotfirstlikelihood(x, x̂, chain; filename="extrapolations", plotpath="", xs=-1.:0.005:1, ys=-1.0:0.005:1.0)
    
    likelihood = (x₁, x₂) -> pdf(MvNormal(x̂[1, :], mean(chain, :Dₓ)), [x₁, x₂])

    xlims = extrema(x̂[:, 1])
    ylims = extrema(x̂[:, 2])
    
    xs = range(xlims..., length=1_000)
    ys = range(ylims..., length=1_000)
    
    heatmap(
        xs, ys, likelihood, 
        title="Position extrapolation with first point likelihood",
        xlabel="x", ylabel="y", legend=:none, aspect_ratio=1,
        xlims=xlims, ylims=ylims, dpi=300)

    scatter!(
        x[:, 1], x[:, 2], 
        label="estimate", legend=:bottomright,
        c=:white, markersize=2)

    scatter!(
        [x̂[1, 1]], [x̂[1, 2]],
        label="extrapolation", c=:black,
        markersize=2)

    savefig(joinpath(plotpath, filename))

end

function plotmontecarlo(x, x̂; plotpath="", filename="mcextrapolation")

    Nsimulations = size(x̂, 3)

    xlims = extrema(vcat(x[:, 1], vec(x̂[1, :, :])))
    ylims = extrema(vcat(x[:, 2], vec(x̂[2, :, :])))

    plot(
        x[:, 1], x[:, 2], label="estimate", c=:blue,
        xlabel="x", ylabel="y", legend=:bottomright, aspect_ratio=1, dpi=500, 
        xlims=xlims, ylims=ylims
    )

    
    for n in 1:Nsimulations

        plot!(
            x̂[1, :, n], x̂[2, :, n], label=:none,
            c=:red, alpha=0.05, aspect_ratio=1,
            xlims=xlims, ylims=ylims
            )
            
    end
    
    x̂ₘ = squeeze(mean(x̂, dims=3))'
        
    plot!(x̂ₘ[:, 1], x̂ₘ[:, 2], label="mean extrapolation", c=:black)
        
    savefig(joinpath(plotpath, filename))

end
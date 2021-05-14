function plotfirstlikelihood(x, x̂, chain; plotpath="", xs=-1.:0.005:1, ys=-1.0:0.005:1.0)
    
    likelihood = (x₁, x₂) -> pdf(MvNormal(x̂[1, :], mean(chain, :Dₓ)), [x₁, x₂])
    
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

function plotmontecarlo(x, x̂; plotpath="")

    Nsimulations = size(x̂, 3)

    lims = extrema(x) .* 1.05

    plot(
        x[:, 1], x[:, 2], label="estimate", c=:blue,
        xlabel="x", ylabel="y", legend=:bottomright, aspect_ratio=1, dpi=300, xlims=lims, ylims=lims
    )


    for n in 1:Nsimulations

        plot!(
            x̂[1, :, n], x̂[2, :, n], label=:none,
            c=:red, alpha=0.01, aspect_ratio=1,
            xlims=lims, ylims=lims
        )

    end

    savefig(joinpath(plotpath, "mcextrapolation"))

end
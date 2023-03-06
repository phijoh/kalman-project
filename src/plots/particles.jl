function plotparticledensity(p::Matrix{Int64}, w::Vector{Float64}, fr::Frame; plotargs...)
    figure = plot()
    plotparticledensity!(figure, p, w, fr; plotargs...)

    return figure
end
function plotparticledensity!(fig, p::Matrix{Int64}, w::Vector{Float64}, fr::Frame; plotargs...)

    width, height = size(fr)
    N = size(p, 1)

    y, x, u, v = getexpectedposition(p, w)
    Σ = getvariance(p, w)

    N = MvNormal([x, y], Σ[1:2, 1:2])
    xs = 1:1:width
    ys = 1:1:height
    f̄ = pdf(N, [x, y])

    heatmap!(fig,
        xs, ys, 1 .- fr;
        c = :grays,
        aspect_ratio=1,
        xlabel="\$x\$", ylabel="\$y\$",
        xlims = extrema(xs), ylims = extrema(ys),
        plotargs...
    )

    contour!(fig,
        xs, ys, (x, y) -> pdf(N, [x, y]);
        alpha = 1., clims = (0, f̄)
    )

    return fig
end

function scatterparticles(particles, fr; weights=nothing, kwargs...)
    N = size(particles, 1)
    S = size(fr, 1)

    alphas = isnothing(weights) ? ones(N) : weights ./ 2maximum(weights)

    fig = heatmap(fr; aspect_ratio=1, xlims=(0, S), ylims=(0, S), color=:greys, kwargs...)

    scatter!(fig, particles[:, 2], particles[:, 1], label=nothing, color=:darkred, alpha=alphas)

    return fig
end


function plotprecision(
    results, duration;
    labels=["speed", "opacity", "noise", "delay"],
    kwargs...
)

    S, T = size(results[:particles])[1:2]
    idxlabels = checkdifferencelabel(results[:specs])

    σ = Matrix{Float64}(undef, T, S) # FIXME: We only use after duration anyway

    for s in 1:S

        particlesovertime = @view results[:particles][s, :, :, :]
        weightsovertime = @view results[:weights][s, :, :]

        for t in 1:T
            p = particlesovertime[t, :, 1:2]
            w = StatsBase.weights(weightsovertime[t, :])

            σ[t, s] = 1 / (var(p[:, 1], w) + var(p[:, 2], w))
        end
    end

    varfig = plot(;
        xlabel="\$t\$",
        ylabel="\$(\\sigma_{x}^2 + \\sigma_{y}^2)^{-1}\$",
        legend=:top,
        kwargs...
    )

    plott = 1:T

    for s in 1:S

        precision = σ[plott, s]

        label = join(
            ("$(t[1]) = $(t[2])"
             for t ∈ zip(labels[idxlabels], results[:specs][s][idxlabels])),
            ", "
        )

        plot!(varfig, plott, precision; label=label)

    end

    vline!(varfig, [duration - 1], linecolor=:black, linestyle=:dot, label=nothing)


    return varfig

end

function plotposition(
    results, duration, dim::Int64;
    labels=["speed", "opacity", "noise", "delay"],
    kwargs...
)
    S, T = size(results[:particles])[1:2]
    idxlabels = checkdifferencelabel(results[:specs])

    X = Matrix{Float64}(undef, T, S) 
    σ = Matrix{Float64}(undef, T, S)

    for s in 1:S
        particlesovertime = @view results[:particles][s, :, :, :]
        weightsovertime = @view results[:weights][s, :, :]

        for t in 1:T
            p = particlesovertime[t, :, :]
            w = weightsovertime[t, :]

            X[t, s] = getexpectedposition(p, w)[dim]
            σ[t, s] = var(p[:, dim], StatsBase.Weights(w)) / sqrt.(length(p[:, dim]))
        end
    end

    varfig = plot(;xlabel="\$t\$", ylabel = "Estimation error", legend=:top, kwargs...)
    plott = 1:T

    for s in 1:S
        speed = results[:specs][s][1]
        
        wedgeposition = (500 - 21) .- range(0; step = speed, length = T)

        label = join(
            ("$(t[1]) = $(t[2])"
            for t ∈ zip(labels[idxlabels], results[:specs][s][idxlabels])),
            ", "
        )

        plot!(varfig, plott, X[:, s]; label=label, ribbon = σ[:, s], fillalpha = 0.1)

    end

    vline!(varfig, [duration - 1], linecolor=:black, linestyle=:dot, label=nothing)

    return varfig

end
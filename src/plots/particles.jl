function plotparticledensity(particles::Matrix{Int64}, weights::Vector{Float64}, fr::Frame; rfsize=2)

    width = size(fr, 1)
    N = size(particles, 1)

    particle_density = zeros(Float64, size(fr))

    for i in 1:N
        x, y = particles[i, 1:2]

        x₀, y₀ = max.([x, y] .- (rfsize - 1), 1)
        x₁, y₁ = min.([x, y] .+ (rfsize - 1), width)

        particle_density[x₀:x₁, y₀:y₁] .+= weights[i]
    end

    figure = contourf(
        particle_density;
        size=(400, 400), legend=nothing, aspect_ratio=1,
        xlabel=L"x", ylabel=L"y", c=:Blues
    )

    return figure

end

function scatterparticles(particles, fr; weights=nothing, kwargs...)
    N = size(particles, 1)
    S = size(fr, 1)

    alphas = isnothing(weights) ? ones(N) : weights ./ 2maximum(weights)

    fig = heatmap(fr; aspect_ratio=1, xlims=(0, S), ylims=(0, S), color=:greys, kwargs...)

    scatter!(fig, particles[:, 2], particles[:, 1], label=nothing, color=:darkred, alpha=alphas)

    return fig
end

function plotexpectedposition(particles, weights, fr; kwargs...)
    S = size(fr, 1)
    x, y = getexpectedposition(particles, weights)
    Σ = sqrt.(getvariance(particles, weights))

    N = MvNormal([x, y], Σ[1:2, 1:2])
    pdfN(x, y) = pdf(N, [x, y])

    xs = ys = 0:1:S

    fig = heatmap(
        xs, ys, fr';
        c=:grays, aspect_ratio=1,
        colorbar=false,
        xlims=ylims = (1, S), kwargs...
    )

    scatter!(
        fig,
        [x], [y];
        markersize=2, markerstrokewidth=0,
        label=nothing, c=:darkred
    )

    contour!(
        fig,
        xs, ys, (a, b) -> pdfN(a, b) / pdfN(x, y);
        clim=(0, 1), c=:Reds,
        colorbar=false
    )

    return fig
end

function plotangle(
    results, duration;
    labels=["speed", "opacity", "noise", "delay"],
    kwargs...)

    S, T = size(results[:particles])[1:2]
    idxlabels = checkdifferencelabel(results[:specs])

    θ = Matrix{Float64}(undef, T, S)


    for s in 1:S

        particlesovertime = @view results[:particles][s, :, :, :]
        weightsovertime = @view results[:weights][s, :, :]

        for t ∈ 1:T
            p = particlesovertime[t, :, 1:2]
            θₜ = xytoangle.(eachrow(p))
            w = weightsovertime[t, :]

            θ[t, s] = θₜ'w
        end
    end

    varfig = plot(
        xlabel=L"t", ylabel=L"\theta_T - \theta_t";
        legend=:bottomleft, kwargs...)


    for s in 1:S

        label = join(
            ("$(t[1]) = $(t[2])"
             for t ∈ zip(labels[idxlabels], results[:specs][s][idxlabels])),
            ", "
        )


        plot!(varfig, 1:T, θ[:, s]; label=label, kwargs...)
    end

    vline!(varfig, [duration - 1], linecolor=:black, linestyle=:dot, label=nothing)

    return varfig
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
            # θₜ = xytoangle.(eachrow(p)) .* (180 / π)
            w = StatsBase.weights(weightsovertime[t, :])

            σ[t, s] = 1 / (var(p[:, 1], w) + var(p[:, 2], w))
        end
    end

    varfig = plot(;
        xlabel=L"t",
        ylabel=L"(\sigma_{x}^2 + \sigma_{y}^2)^{-1}",
        legend=:top,
        kwargs...
    )

    plott = 1:T

    for s in 1:S

        precision = 1 ./ σ[plott, s]

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
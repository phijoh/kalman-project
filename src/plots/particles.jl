function plotparticledensity(particles::Matrix{Int64}, weights::Vector{Float64}, frame::Frame; plotargs...)
    fig = plot()
    plotparticledensity!(fig, particles, weights, frame; plotargs...)

    return fig
end
function plotparticledensity!(fig, particles::Matrix{Int64}, weights::Vector{Float64}, fr::Frame; plotargs...)

    width, height = size(fr)
    N = size(particles, 1)

    x, y, u, v = getexpectedposition(particles, weights)
    Σ = getvariance(particles, weights)

    N = MvNormal([x, y], Σ[1:2, 1:2])
    xs = 1:1:width
    ys = 1:1:height
    f̄ = pdf(N, [x, y])

    plotframe!(fig, fr;
        c = :grays,
        aspect_ratio=1,
        xlabel="\$x\$", ylabel="\$y\$",
        xlims = extrema(xs), ylims = extrema(ys),
        plotargs...
    )

    contour!(fig,
        xs, ys, (x, y) -> pdf(N, [x, y]) / f̄;
        alpha = 1.
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

function plotposition(
    results, models::Vector{TwinkleGoesParameters}; 
    variableindex = 1, # x = 1 and y = 2
    labels, kwargs...)
    
    model = first(models)
    Tᵢ = mstoframes(model.inducerduration)
    Tₙ = mstoframes(model.noiseduration)
    T = Tᵢ + Tₙ

    S = length(models)

    labels = isnothing(labels) ? ["model $i" for i in 1:S] : labels

    X = Matrix{Float64}(undef, T, S) 
    σ = Matrix{Float64}(undef, T, S)

    for s ∈ 1:S
        particlesovertime, weightsovertime, frames = results[s]
        model = models[s]


        for t in 1:T
            p = particlesovertime[t, :, :]
            w = weightsovertime[t, :]

            X[t, s] = getexpectedposition(p, w)[variableindex]
            σ[t, s] = var(p[:, variableindex], StatsBase.Weights(w)) / sqrt.(length(p[:, variableindex]))
        end
    end

    wedgeposition = Vector{Float64}(undef, T)
    wedgeposition[1:Tᵢ] = range(21; step = first(models).speed, length = Tᵢ)
    wedgeposition[(Tᵢ + 1):end] .= wedgeposition[Tᵢ]

    plott = framestoms.(1:T)

    varfig = plot(plott, wedgeposition; xlabel="\$t\$", ylabel = "Estimated position", legend=:topleft, label = "Wedge position", c = :black, linestyle = :dash, ylims = (0, Inf), kwargs...)

    for s in 1:S        
        plot!(varfig, plott, X[:, s]; label=labels[s], ribbon = σ[:, s], fillalpha = 0.1)
    end

    vline!(varfig, [framestoms(Tᵢ)], linecolor=:black, linestyle=:dot, label=nothing)

    return varfig

end

function plotvelocity(
    results, models::Vector{TwinkleGoesParameters}; 
    variableindex = 1, # x = 1 and y = 2
    T₀ = 200,
    labels, kwargs...)
    
    model = first(models)
    Tᵢ = mstoframes(model.inducerduration)
    Tₙ = mstoframes(model.noiseduration)
    T = Tᵢ + Tₙ

    S = length(models)

    labels = isnothing(labels) ? ["model $i" for i in 1:S] : labels

    X = Matrix{Float64}(undef, T, S) 
    σ = Matrix{Float64}(undef, T, S)

    for s ∈ 1:S
        particlesovertime, weightsovertime, frames = results[s]
        model = models[s]


        for t in 1:T
            p = particlesovertime[t, :, :]
            w = weightsovertime[t, :]

            X[t, s] = getexpectedposition(p, w)[2 + variableindex]
            σ[t, s] = var(p[:, variableindex], StatsBase.Weights(w)) / sqrt.(length(p[:, variableindex]))
        end
    end

    wedgevelocity = Vector{Float64}(undef, T)
    wedgevelocity[1:Tᵢ] .= speed
    wedgevelocity[(Tᵢ + 1):end] .= 0.

    plott = framestoms.(1:T)

    varfig = plot(plott, wedgevelocity; xlabel="\$t\$", ylabel = "Estimated position", legend=:topleft, label = "Wedge velocity", c = :black, linestyle = :dash, kwargs...)

    for s in 1:S        
        plot!(varfig, plott, X[:, s]; label=labels[s], ribbon = σ[:, s], fillalpha = 0.1)
    end

    vline!(varfig, [framestoms(Tᵢ)], linecolor=:black, linestyle=:dot, label=nothing)

    return varfig

end
function plotparticledensity(particles::Matrix{Int64}, weights::Vector{Float64}, fr::Frame; rfsize = 2)

    width = size(fr, 1)
    N = size(particles, 1)

    particle_density = zeros(Float64, size(fr))

    for i in 1:N
        x, y = particles[i, 1:2]

        x₀,y₀ = max.([x, y] .- (rfsize-1), 1)
        x₁,y₁ = min.([x, y] .+ (rfsize-1), width)

        particle_density[x₀:x₁, y₀:y₁] .+= weights[i]
    end

    figure = contourf(
        particle_density; 
        size = (400, 400), legend = nothing, aspect_ratio = 1,
        xlabel = L"x", ylabel = L"y", c = :Blues
    )

    return figure

end

function scatterparticles(particles, fr)
    S = size(fr, 1)
    fig = heatmap(fr; aspect_ratio = 1, xlims = (0, S), ylims = (0, S), color = :greys)

    scatter!(fig, particles[:, 2], particles[:, 1], label = nothing, color = :white)

    return fig
end

function plotexpectedposition(particles, weights, fr)
    S = size(fr, 1)
    fig = heatmap(fr; aspect_ratio = 1, xlims = (0, S), ylims = (0, S), color = :greys)

    x = mean(particles[:, 1], StatsBase.weights(weights))
    y = mean(particles[:, 2], StatsBase.weights(weights))

    scatter!(fig, [x], [y], label = nothing, color = :red)

    return fig
end

function plotvariance(particlesovertime, weightsovertime, duration; kwargs...)
    T = size(weightsovertime, 1)
    σ = Matrix{Float64}(undef, T, 2)

    for t in 1:T
        p = particlesovertime[t, :, 1:2]
        w = StatsBase.weights(weightsovertime[t, :])

        σ[t, 1] = var(p[:, 1], w)
        σ[t, 2] = var(p[:, 2], w)
    end

    varfig = plot(xlabel = L"t", ylabel = L"\sigma^2"; legend = :bottomleft, kwargs...)
    plot!(varfig, 1:T, σ[:, 1]; label = L"\sigma^2_x")
    plot!(varfig, 1:T, σ[:, 2]; label = L"\sigma^2_y")
    vline!(varfig, [duration]; linestyle = :dash, c = :black, label = "Disappearance")


    return varfig
end
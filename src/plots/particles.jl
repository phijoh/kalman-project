function plotparticles(particles::Matrix{Int64}, fr::Matrix{Float64})

    particle_density = zeros(Int64, size(fr))

    for (x, y, _, _) in eachrow(particles)
        particle_density[x, y] += 1
    end

    figure = heatmap(particle_density, size = (400, 400), legend = nothing, aspect_ratio = 1)

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
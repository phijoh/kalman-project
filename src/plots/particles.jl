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

function getexpectedvalue(p, w)

    x = mean(p[:, 1], StatsBase.weights(w))
    y = mean(p[:, 2], StatsBase.weights(w))

    θ = xytoangle([x y]) .* (180 / π)

    return θ
end


function plotangle(results, duration; kwargs...)
    S, T, _, _ = size(results[:particles])

    σ = Matrix{Float64}(undef, T, S) # FIXME: We only use after duration anyway
    θ = Matrix{Float64}(undef, T, S)

    for s in 1:S

        particlesovertime = @view results[:particles][s, :, :, :]
        weightsovertime = @view results[:weights][s, :, :]

        θ₀ = getexpectedvalue(
            results[:particles][s, duration, :, :],
            results[:weights][s, duration, :]
        )

        for t in 1:T
            p = particlesovertime[t, :, 1:2]
            θₜ = xytoangle.(eachrow(p)) .* (180 / π)
            w = StatsBase.weights(weightsovertime[t, :])

            θ[t, s] = θ₀ - θₜ'w
            σ[t, s] = std(θₜ, w) 
        end
    end

    varfig = plot(
        xlabel = L"t", ylabel = L"\theta"; 
        legend = :bottomleft, kwargs...)

    plott = 1:T
    
    for s in 1:S
        plot!(
            varfig, plott, θ[plott, s];
            ribbon = σ[plott, s], fillalpha = 0.2,
            label = "speed = $(results[:speeds][s])")
    end

    vline!(varfig, [duration], linecolor = :black, linestyle = :dot,
    label = "wedge duration")    
    # vline!(varfig, [duration]; linestyle = :dash, c = :black, label = "Disappearance")

    return varfig
end



function plotprecision(results, duration; kwargs...)

    S, T, _, _ = size(results[:particles])

    σ = Matrix{Float64}(undef, T, S) # FIXME: We only use after duration anyway

    for s in 1:S

        particlesovertime = @view results[:particles][s, :, :, :]
        weightsovertime = @view results[:weights][s, :, :]

        for t in 1:T
            p = particlesovertime[t, :, 1:2]
            θₜ = xytoangle.(eachrow(p)) .* (180 / π)
            w = StatsBase.weights(weightsovertime[t, :])

            σ[t, s] = var(p[:,1], w) + var(p[:,2], w)
        end
    end

    varfig = plot(
        xlabel = L"t", ylabel = L"(\sigma_{x}^2 + \sigma_{y}^2)^{-1}"; 
        legend = :outerright, kwargs...)

    plott = 1:T
    
    for s in 1:S
        plot!(
            varfig, plott, 1 ./σ[plott, s];
            label = "speed = $results[:specs][s][1], 
            opacity =  $results[:specs][s][2],
            dynamic = $results[:specs][s][3]")
    end

    vline!(varfig, [duration-1], linecolor = :black, linestyle = :dot,
    label = "stimulus duration")    
    # vline!(varfig, [duration]; linestyle = :dash, c = :black, label = "Disappearance")

    return varfig

end
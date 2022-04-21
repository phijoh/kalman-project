xytoangle(coord) = atan(coord[2] / coord[1])

function getexpectedposition(p::Matrix{Int64}, w::Vector{Float64})

    x = mean(p[:, 1], StatsBase.weights(w))
    y = mean(p[:, 2], StatsBase.weights(w))

    return x, y
end

function getvariance(p::Matrix{Int64}, w::Vector{Float64})
    dimensions = size(p, 2)
    Σ = zeros(dimensions, dimensions)

    for d ∈ 1:dimensions
        Σ[d, d] = var(p[:, d], StatsBase.weights(w))
    end

    return Σ
end

function getexpectedangle(p::Matrix{Int64}, w::Vector{Float64})

    x, y = getexpectedposition(p, w)
    θ = xytoangle([x y])
    return θ
end

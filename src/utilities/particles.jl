function getexpectedposition(p::Matrix{Int64}, w::Vector{Float64})

    x = mean(p[:, 1], StatsBase.weights(w))
    y = mean(p[:, 2], StatsBase.weights(w))
    u = mean(p[:, 3], StatsBase.weights(w))
    v = mean(p[:, 4], StatsBase.weights(w))

    return x, y, u, v
end

function getvariance(p::Matrix{Int64}, w::Vector{Float64})
    dimensions = size(p, 2)
    Σ = zeros(dimensions, dimensions)
    weights = StatsBase.weights(w)

    for d ∈ 1:dimensions
        Σ[d, d] = var(p[:, d], weights)
    end

    return Σ
end

function integernormal(N::Int64, upper::Int64)
    dist = Binomial(upper * 2, 0.5)

    return (rand(dist, N) .- upper) .÷ 2
end
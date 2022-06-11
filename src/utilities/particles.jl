function computeangle(x, y)
    θ = atan(y, x)
    θ < 0 ? π + abs(θ) : θ
end
xytoangle(coord) = computeangle(coord[1], coord[2])
uvtoangle(coord) = computeangle(coord[3], coord[4])

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
    getexpectedposition(p, w) |> xytoangle
end


function getcoherencevelocity(particles, weights)
    angles = map(uvtoangle, eachrow(particles))
    return mean(angles, StatsBase.weights(weights))
end

function integernormal(N::Int64, upper::Int64)
    dist = Binomial(upper * 2, 0.5)

    return rand(dist, N) .- upper
end
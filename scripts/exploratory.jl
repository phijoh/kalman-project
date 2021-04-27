using Distributions, Random

using StateSpaceRoutines

Random.seed!(121212)

σy = √(π)
τ = π/5
T = 100

ε = Normal(0, σy^2)
Y = range(0, τ*T, length = T) .+ (rand(ε, T) .* range(1.0, 0.01, length=T))

P₁ = 3.5
a₁ = 0

A = zeros(Float64, T)
P = copy(A)
A[1] = a₁
P[1] = P₁

for (t, y) in enumerate(Y)
    if t == T continue end

    a, p = A[t], P[t]

    v = y - a # Error
    F = σy^2 + p # Variance update

    K = p / F

    p′ = p*(1-K) + 1
    a′ = a + K*v

    A[t+1] = a′
    P[t+1] = p′ + 1

    gain[t] = K

    print("Step $t: v = $v, τ ≈ $K \n")

end

plot(1:T, A, ribbon = P, label= "filter", legend=:bottomright)
scatter!(1:T, Y, label = "data", markersize=1.2)
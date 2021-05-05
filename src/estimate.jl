Random.seed!(12)

@model movement(x, v, v₁, Δ) = begin

    T, M = size(v)

    k = 1.
    aₓ ~ InverseGamma(2, 3)
    h ~ InverseGamma(2, 3)

    σₓ = k + aₓ
    σᵧ = h * σₓ

    Σ = [
        σₓ 0 ;
        0 σᵧ
    ]

    v[1, :] ~ MvNormal(v₁, Σ)

    for t = 2:T v[t, :] ~ MvNormal(v[t - 1, :], Σ) end

end

"""
Generates data with T time points and two variances
"""
function dgp(σₓ, σᵧ; T=200)
    
    Σ = [σₓ 0; 0 σᵧ]

    vel = zeros(T, 2)

    for t in 2:T
        veldist = MvNormal(vel[t - 1, :], Σ)
        vel[t, :] = rand(veldist)
    end

    return vel

end
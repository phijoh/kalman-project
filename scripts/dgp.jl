"""
Generates data with pixels and T time points given a multivariate normal ε
"""
function generatedata(pixels::Int64, T::Int64, ε::MultivariateNormal; dθ=π / 64)

    radius = pixels ÷ 5
    θs = range(0., length=T, step=dθ)

    noise = rand(ε, T)'

    xs = @. radius * sin(θs) 
    ys = @. radius * cos(θs) 

    Y = zeros(Int64, T, 2)
    Y[:, 1] =  round.(Int64, xs)
    Y[:, 2] = round.(Int64, ys)

    return Y .+ round.(Int64, noise)

end
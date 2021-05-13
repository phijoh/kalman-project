Id(M) = Matrix{Float64}(I, M, M)

"""
Drop singleton dimensions
"""
function squeeze(A) 
    dropdims(A, dims=(findall(size(A) .== 1)...,))
end
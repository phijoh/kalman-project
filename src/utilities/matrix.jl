"""
Drop singleton dimensions
"""
function squeeze(A) 
    dropdims(A, dims=(findall(size(A) .== 1)...,))
end


function normalise(x)
    x ./ sum(x)
end
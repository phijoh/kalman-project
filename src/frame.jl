struct Frame <: AbstractMatrix{Float64}
    M::Matrix{Float64}
end


Base.size(F::Frame, dims...) = size(F.M, dims...)

function coordtoindex(F::Frame, I::Vararg{Int,2})
    x, y = I
    column = x
    row = size(F, 1) - (y - 1)

    return row, column
end

function Base.getindex(F::Frame, I::Vararg{Int,2})
    row, column = coordtoindex(F, I...)
    getindex(F.M, row, column)
end
function Base.setindex!(F::Frame, v, I::Vararg{Int,2})
    row, column = coordtoindex(F, I)
    F.M[row, column] = v
end
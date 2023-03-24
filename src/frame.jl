"""
Subtype of matrix, that can be indexed with (x, y) and returns (i, j) = (n - y, x). This allows us to work with coordinates in the frame insteady of entries of the matrix.
"""
struct Frame <: AbstractMatrix{Float64}
    M::Matrix{Float64}
end


Base.size(F::Frame, dims...) = size(F.M', dims...)

coordtoindex(F::Frame, I::Tuple{Int64, Int64}) = coordtoindex(F, I...)
function coordtoindex(F::Frame, I::Vararg{Int,2})
    x, y = I
    column = x
    row = size(F, 2) - (y - 1)

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

"""
Redefine product between to frames to be element wise
"""
function Base.:*(A::Frame, B::Frame)
    Frame(A.M .* B.M)
end

function Base.:+(A::Frame, B::Frame)
    Frame(A.M + B.M)
end
function Base.:-(A::Frame, B::Frame)
    Frame(A.M - B.M)
end

function Base.:*(s::Number, F::Frame)
    Frame(s * F.M)
end

function Base.copy(F::Frame)
    Frame(copy(F.M))
end

function translateframe(frame::Frame, Δx::Int64, Δy::Int64)
    Frame(warp(frame.M, Translation(-Δy, -Δx), indices_spatial(frame.M), 0))
end

function randomframe(widthframe, heightframe)::Frame
    Frame(rand(heightframe, widthframe))
end
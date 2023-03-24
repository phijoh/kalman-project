function mstoframes(ms::Real)::Int64
    ceil(Int64, ms * framespersecond / 1000)
end

function framestoms(frames::Int64)
    frames * (1000 ÷ framespersecond)
end
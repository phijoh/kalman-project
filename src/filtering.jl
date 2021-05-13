
μ₀ = zeros(2)

# FIXME
function cartesiantomatrix(indices::Vector{CartesianIndex{2}})
    N = length(indices)
    mform = zeros(Int64, N, 2)

    for (t, coord) in enumerate(indices)
        x, y = Tuple(coord)
        mform[t, :] = [x, y]
    end

    return mform
end

function detect(frame::Matrix{Float64})
    imgg = imfilter(frame, Kernel.LoG(25))

    
    lb = minimum(imgg)
    dark = findall(c -> c ≈ lb, imgg) |> cartesiantomatrix
    y, x = mean(dark, dims=1)

    return x, y
end

function filtering(frames)

    Nframes, height, width = size(frames)

    position = zeros(Nframes, 2)

    @threads for frame in 1:Nframes
        x, y = detect(frames[frame, :, :])
        position[frame, :] = [x, y]
    end

    return position

end

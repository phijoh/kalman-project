
μ₀ = zeros(2)

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

    Nframes  = size(frames, 1)
    side = size(frames, 2)

    position = zeros(Nframes, 2)

    @threads for f in 1:Nframes
        x, y = detect(frames[f, :, :])
        position[f, :] = @. (2 * [x, y] / side) - 1
    end

    return position

end

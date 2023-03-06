"""
ms to number of frames
"""
function mstoframes(ms)
    ms * framespersecond / 1000
end

function framestoms(frames::Int64)
    frames * (1000 ÷ framespersecond)
end

function framegenerator(inducerduration::Int64, noiseduration::Int64, speed::Int64, opacity::Float64; dynamic=false, framesize=500::Int64, wedgesize=42::Int64)::Vector{Frame}

    wedge = zeros(Bool,framesize,framesize)
    wedge[1:wedgesize,round(Int,framesize/2-wedgesize/2):round(Int,framesize/2+wedgesize/2)] .= 1
    
    currentwedge = opacity * copy(wedge)
    inducernoise = rand(framesize, framesize)
    alphablend(frame) = @. frame + inducernoise * (1 - frame) 

    inducerframes = ceil(Int64, mstoframes(inducerduration))
    noiseframes = ceil(Int64, mstoframes(noiseduration))

    T = inducerframes + noiseframes
    frames = zeros(T, framesize, framesize)
    frames[1, :, :] = alphablend(currentwedge)

    t = Translation(-speed, 0)
    
    for inducerframe in 2:inducerframes 

        currentwedge = warp(currentwedge, t, indices_spatial(currentwedge), 0)

        frames[inducerframe, :, :] = alphablend(currentwedge)

    end

    @threads for noiseframe in inducerframes:T 
        frames[noiseframe, :, :] = dynamic ? 
            rand(framesize, framesize) : 
            inducernoise 
    end

    return [Frame(frames[t, :, :]) for t ∈ 1:T]

end 
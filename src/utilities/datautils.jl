rgbtogrey(frames) = mean(frames; dims=4)
scale(frames) = frames ./ 255

"""
Reads the stimulus .mat file
"""
function ingeststimulus(filepath)
    vars = matread(filepath)

    if "frames" in keys(vars)

        greyframes = vars["frames"] |> rgbtogrey |> squeeze |> scale
        showframe = greyframes[2:end, :, :] # Skip title frame

        return showframe

    else
        throw("No frames found in the .mat file")
    end
end

"""
Compute velocity given position matrix and time step
"""
function getvelocity(x::Matrix{Float64}, Δt::Float64)
    T, M = size(x)
    u = zeros(T, M)

    @threads for t in 1:(T - 1) 
        u[t, :] = (x[t + 1, :] - x[t, :]) / Δt 
    end

    return u
end

"""
Make constructor for frames from TG_test.mat
"""
function framegeneratorfactory(tgfile)
    
    wedge = tgfile["stimMask"]
    flashduration = convert(Int64, tgfile["flashDur"])

    S = size(wedge, 1)

    function framegenerator(speed::Float64, inducerduration::Int64, opacity::Float64; dynamic=false)::Vector{Frame}

        Δθ = deg2rad(speed)
        
        currentwedge = opacity * copy(wedge)
        inducernoise = rand(S, S)
        alphablend(frame) = @. frame + inducernoise * (1 - frame) 

        T = inducerduration + flashduration
        frames = zeros(T, S, S)
        frames[1, :, :] = alphablend(currentwedge)
        
        for inducerframe in 2:inducerduration 

            currentwedge = imrotate(currentwedge, Δθ, axes(currentwedge))
            replace!(currentwedge, NaN => 0.)

            frames[inducerframe, :, :] = alphablend(currentwedge)

        end

        @threads for noiseframe in inducerduration:T frames[noiseframe, :, :] = dynamic ? rand(S, S) : inducernoise end

        return [Frame(frames[t, :, :]) for t ∈ 1:T]

    end 

    return framegenerator

end

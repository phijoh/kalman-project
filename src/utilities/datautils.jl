rgbtogrey(frames) = mean(frames; dims=4)
scale(frames) = frames ./ 255

"""
Reads the stimulus .mat file
"""
function ingeststimulus(filepath)
    vars = matread(filepath)

    if "frames" in keys(vars)

        frames = vars["frames"]         
        greyframes = frames |> rgbtogrey |> squeeze |> scale


        showframe = greyframes[2:end, :, :] # Skip title frame

        return showframe

    else

        throw("No frames found in the .mat file")
         
    end
end

function getvelocity(x::Matrix{Float64})
    T, M = size(x)
    u = zeros(T, M)

    for t in 2:T u[t, :] = x[t, :] - x[t - 1, :] / Δt end

    return u
end
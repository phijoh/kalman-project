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
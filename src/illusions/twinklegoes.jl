struct TwinkleGoesParameters
    speed::Int64 # in frames per ms
    opacity::Real
    noise::Bool
    inducerduration::Real # in ms
    noiseduration::Real # in ms
end

function generateframes(
    model::TwinkleGoesParameters; 
    widthframe = 500, heightframe = 500, wedgesize = 42)::Vector{Frame}

    wedge = Frame(zeros(heightframe, widthframe))
    wedge[
        1:wedgesize, # x coordinates
        (heightframe - wedgesize) ÷ 2:(heightframe + wedgesize) ÷ 2 # y coordinates
    ] .= model.opacity
        
    inducernoise = randomframe(widthframe, heightframe) # noise in the inducer phase

    Tᵢ = mstoframes(model.inducerduration)
    Tₙ = mstoframes(model.noiseduration)

    currentwedge = copy(wedge)
    frames = [wedge * inducernoise]

    for t in 2:Tᵢ 
        currentwedge = translateframe(currentwedge, model.speed, 0)
        push!(frames, currentwedge * inducernoise)
    end

    for t in Tᵢ:(Tᵢ + Tₙ) 
        push!(frames, model.noise ? randomframe(widthframe, heightframe) : inducernoise)
    end

    return frames 
end
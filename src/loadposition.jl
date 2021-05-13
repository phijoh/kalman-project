function loadposition(datapath; cache=false, verbose=false, framelimit=nothing, framesize=750)

    inpath = joinpath(datapath, "stimulus_frames.mat")
    outpath = joinpath(datapath, "position.jld")

    if cache && isfile(outpath)
        verbose && println("Using cached position file $outpath...")

        x = load(outpath, "position")
    else
        verbose && println("Loading frames from $inpath...")
        frames = ingeststimulus(inpath)

        T, width, height = size(frames)

        hcenter = width ÷ 2
        vcenter = height ÷ 2

        hsize = vsize = framesize ÷ 2 # FIXME: Account for non square screen

        hl, hu = hcenter - hsize, hcenter + hsize - 1
        vl, vu = vcenter - vsize, vcenter + vsize - 1

        upper = isnothing(framelimit) ?  T : framelimit
    
        activeframes = frames[2:upper, hl:hu, vl:vu]

        verbose && println("Computing wedge position...")
    
        x = filtering(activeframes)

        save(outpath, "position", x)
    end


    return x

end
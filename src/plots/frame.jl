"""
Plot frame assuming top-left coordinate system.
"""
function plotframe(frame::Frame; plotkwargs...)
    fig = plot(); plotframe!(fig, frame; plotkwargs...)
    return fig
end
function plotframe!(fig, frame::Frame; plotkwargs...)
    heatmap!(fig, 1:size(frame, 1), 1:size(frame, 2), frame'; yflip = true, c = :grays, colorbar = false, plotkwargs...)
end
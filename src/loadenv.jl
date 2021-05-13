using DotEnv

# Environment variables
DotEnv.config()
shallplot = getbool(ENV, "PLOT")
verbose = getbool(ENV, "VERBOSE")
cache = getbool(ENV, "CACHE")

framespersecond = getint(ENV, "FRAMES_PER_SECOND"; def="60")
seed = getint(ENV, "RANDOM_SEED")

datapath = get(ENV, "INPUT", "data")
plotpath = get(ENV, "PLOT_PATH", nothing)
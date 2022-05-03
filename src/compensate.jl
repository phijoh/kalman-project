function compensate(particles, τ, framesize)

    compensatedparticles = copy(particles)
    dims = size(particles,2)
    for _ in 1:τ
        moveparticles!(compensatedparticles, zeros(dims,dims), framesize)
    end

    return compensatedparticles

end
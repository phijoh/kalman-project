function kalmanfilter(Y; a₁=[0, 0], Σ=Matrix(25 * I, 2, 2))

    T, D = size(Y)

    A = zeros(Int64, T, D)
    A[1, :] = a₁
    
    P = zeros(Int64, T, D, D)
    P[1, :, :] = Matrix(400 * I, D, D)

    posterror = copy(A)

    for t in 1:T - 1

        y = Y[t, :]

        a, p = A[t, :], P[t, :, :] 

        F = p + Σ
        v = y - a # Error

        K = p * inv(F) # Kalman gain

        a′ = round.(Int64, a + K * v)
        p′ = round.(Int64, (I - K) * F)

        A[t + 1, :] = a′
        P[t + 1, :, :] = p′

        posterror[t, :] = y - a′
    end

    return A, P, posterror
end
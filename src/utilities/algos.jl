function checkdifferencelabel(list)

    n = length(first(list))
    seed = [Set() for _ ∈ 1:n]

    for tuple ∈ list
        for (n, elem) ∈ enumerate(tuple)
            push!(seed[n], elem)
        end
    end

    return [length(s) > 1 for s ∈ seed]

end

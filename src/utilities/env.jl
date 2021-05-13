function getbool(env, key)
    (get(ENV, "PLOT", "false") |> lowercase) == "true"
end
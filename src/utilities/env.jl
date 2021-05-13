function getbool(env, key)
    (get(env, key, "false") |> lowercase) == "true"
end

function getint(env, key; def="0")
    strvalue = get(env, key, def)
    return parse(Int64, strvalue) 
end
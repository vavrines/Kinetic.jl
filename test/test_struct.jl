D = read_dict("config.txt")
for key in keys(D)
    s = Symbol(key)
    @eval $s = $(D[key])
end

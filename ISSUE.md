# Known Issues

1. FFT doesn't work well with offset arrays.
1. JLD2 has trouble saving a struct with function fields and the behavior is unpredictable, see https://github.com/JuliaIO/JLD2.jl/issues/401. BSON partially supports this functionality, but the cross-referencing of functions might be a problem.

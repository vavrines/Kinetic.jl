using PyCall

# try importing meshio
# if not existed, install it
try
    meshio = pyimport("meshio")
catch
    using Conda
    Conda.add_channel("conda-forge")
    Conda.add("meshio")
    meshio = pyimport("meshio")
end
include("groups.jl")
include("VecGtensor.jl")
include("initialization.jl")
include("renormalization.jl")
include("physicaldata.jl")


model = "Ising"
Dcut = 16
n = 20
T = initialization(model)
T_fix = RG(T, Dcut, n)
energy = spectrum(T_fix)
using Documenter
include("../src/main.jl")
using .VecG_TNR

makedocs(
    # strict = false, 
    modules = [VecG_TNR],
    sitename = "VecG_TNR documentation",
)
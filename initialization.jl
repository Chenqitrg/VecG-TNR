function initialization(model::AbstractString)
    if model == "Ising"
        Z2 = CyclicGroup(2)
        e = GroupElement(0, Z2)
        a = GroupElement(1, Z2)
        object = Object(Dict(e=>1,a=>1),Z2)
        return random_block_tensor(Float64, Z2, (object,object,object,object))
    end
end
using Revise
includet("main.jl")
using .VecG_TNR

function test_sectoutin()
    D6 = DihedralGroup(3)
    e = identity_element(D6)
    s = GroupElement((1,0), D6)
    r = GroupElement((0,1), D6)
    A = Obj(e=>2, s=>3, r=>2, s*r=>15) 
    B = Obj(e=>2, s*r*r=>4, r=>3, s*r=>2)
    C = Obj(e=>2, s=>3, r*r=>2, s*r=>2)
    D = Obj(e=>2, s=>4, r=>3, s*r*r=>15)
    TT = random_mor(Float64, (A, B, C, D)) # The data is too long to show. Random morphism is generated   
    @show to_sector_outin(keys(TT.data), 2, s)
    # @show to_sector_outin(keys(TT.data), 2, r)
    # @show to_sector_outin(keys(TT.data), 3, s)


end

test_sectoutin()
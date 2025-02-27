using Revise
includet("main.jl")
using .VecG_TNR

function test_group()
    Z2 = CyclicGroup(2)
    @show e = identity_element(Z2) # e
    @show a = GroupElement(1, Z2) # a
    @show a*a # e
    @show inverse(a) # a
    @show a*e == a # true

    Z4 = CyclicGroup(4)
    @show e = identity_element(Z4) # e
    @show a = GroupElement(1, Z4) # a
    @show a*a == GroupElement(2, Z4) # true
    @show inverse(a) # a³

    D6 = DihedralGroup(3) # 𝔻₆
    @show e = identity_element(D6) # e
    @show s = GroupElement((1,0), D6) # s
    @show r = GroupElement((0,1), D6) # r
    @show r*r # r²
    @show r*r*r # r³
    @show r * s # sr²
end

function test_obj()
    D6 = DihedralGroup(3)
    e = identity_element(D6)
    s = GroupElement((1,0), D6)
    r = GroupElement((0,1), D6)
    @show A = Obj(e=>2, s=>3, r=>2, s*r=>15) # 2e ⊕ 3s ⊕ 2r ⊕ 15sr
    @show zero_obj(D6) # nothing
    @show dual_obj(A) # 2e ⊕ 2r² ⊕ 3s ⊕ 15sr
    @show A == dual_obj(dual_obj(A)) # true

end

function test_sector()
    D6 = DihedralGroup(3)
    e = identity_element(D6)
    s = GroupElement((1,0), D6)
    r = GroupElement((0,1), D6)
    @show sect = Sector(e, s, r, s*r)
    # @show Sector(s,s,r,r) should throw an error
    @show sect[1]
    @show sect[2]
    @show sect[3]
    @show sect[4]
end

function test_Mor()
    D6 = DihedralGroup(3)
    e = identity_element(D6)
    s = GroupElement((1,0), D6)
    r = GroupElement((0,1), D6)
    A = Obj(e=>2, s=>3, r=>2, s*r=>15) 
    B = Obj(e=>2, s=>4, r=>3, s*r=>2)
    C = Obj(e=>2, s=>3, r=>2, s*r=>2)
    D = Obj(e=>2, s=>4, r=>3, s*r*r=>15)
    @show T = Mor(Float64, (A, B, C, D))
    @show T.data # Dict{Sector{DihedralGroup}, Array{Float64}}()
    @show T.objects # (2e ⊕ 2r ⊕ 3s ⊕ 15sr, 2e ⊕ 3r ⊕ 4s ⊕ 2sr, 2e ⊕ 2r ⊕ 3s ⊕ 2sr, 2e ⊕ 3r ⊕ 4s ⊕ 15sr²)
    # @show T[e,e,e,e] ERROR: KeyError: key e ⊗ e ⊗ e ⊗ e not found
    @show T[e,e,e,e] = rand(2,2,2,2)
    # @show T[e,e,e,e] = rand(2,3,2,3) ERROR: ArgumentError: Data size (2, 3, 2, 3) does not match sector size (2, 2, 2, 2)
    @show T # The data is too long to show. The e ⊗ e ⊗ e ⊗ e sector is rand(2,2,2,2) given above
    @show T[e,s,s,e]  = rand(2,4,3,2)
    @show T # The data is too long to show. The e ⊗ s ⊗ s ⊗ e sector is rand(2,4,3,2) given above
    @show TT = random_mor(Float64, (A, B, C, D)) # The data is too long to show. Random morphism is generated   
    @show TT = zero_mor(Float64, (A, B, C, D)) # The data is too long to show. Zero morphism is generated
    @show Mat = identity_mor(Float64, A) # The data is too long to show. Identity morphism is generated

    
end

function test_Mor_Z2()
    Z2 = CyclicGroup(2)
    e = identity_element(Z2)
    a = GroupElement(1, Z2)
    B = Obj(e=>1, a=>1)
    C = Obj(e=>1, a=>1)
    @show T = random_mor(Float64, (B,C))
    A = Obj(e=>2)
    @show Mat = identity_mor(Float64, A) # Mor{CyclicGroup, Float64}((2e, 2e), Dict{Sector{CyclicGroup}, Array{Float64}}(e ⊗ e => [1.0 0.0; 0.0 1.0], a ⊗ a => Matrix{Float64}(undef, 0, 0)))
    @show Mat2 = identity_mor(Float64, B) # Mor{CyclicGroup, Float64}((e ⊕ a, e ⊕ a), Dict{Sector{CyclicGroup}, Array{Float64}}(e ⊗ e => [1.0;;], a ⊗ a => [1.0;;]))
end

function test_get_indices()
    D6 = DihedralGroup(3)
    e = identity_element(D6)
    s = GroupElement((1,0), D6)
    r = GroupElement((0,1), D6)
    A = Obj(e=>2, s=>3, r=>2, s*r=>15) 
    B = Obj(e=>2, s=>4, r=>3, s*r=>2)
    C = Obj(e=>2, s=>3, r=>2, s*r=>2)
    D = Obj(e=>2, s=>4, r=>3, s*r*r=>15)
    T = random_mor(Float64, (A, B, C, D))

    @show T[e,e,e,e] # 2×2×2×2 Array{Float64,4}
    @show T[1] # 2e ⊕ 2r ⊕ 3s ⊕ 15sr
    @show T[2] # 2e ⊕ 3r ⊕ 4s ⊕ 2sr
    @show T[2:3] # (2e ⊕ 3r ⊕ 4s ⊕ 2sr, 2e ⊕ 2r ⊕ 3s ⊕ 2sr)
    @show T[end] # 2e ⊕ 3r ⊕ 4s ⊕ 15sr²
    @show lastindex(T) # 4
    @show get_sector_size(T, (e,e,e,e)) # (2, 2, 2, 2)
    @show size = get_sector_size(T, (e,s,s,e)) # (2, 4, 3, 2)
    # @show get_sector_size(T, (e,s,r,e)) ERROR: ArgumentError: The sector (e, s, r, e) is not consistent
    S = Sector(e, s, s, e)
    @show T[S] = rand(size...) # The data is too long to show. Random data is generated and assigned to the sector (e, s, s, e)

end

function test_accend()
    @show is_accend((4,1), 4)
end
# test_group()
# test_obj()
# test_sector()
# test_Mor()
# test_Mor_Z2()
# test_get_indices()
test_accend()
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

    @show Z2xZ4 = ProductGroup(Z2, Z4) # ℤ₂ × ℤ₄
    @show ProductGroup(D6, Z2) # 𝔻₆ × ℤ₂

    aZ2 = GroupElement(1, Z2)
    eZ4 = identity_element(Z4)
    @show ae = GroupElement((aZ2, eZ4), ProductGroup(Z2, Z4)) # (a, e)
    @show sa = GroupElement((s, aZ2), ProductGroup(D6, Z2)) # (s, a)
    @show iter = elements(Z2xZ4)
    @show collect(iter) # Tuple{GroupElement{CyclicGroup}, GroupElement{CyclicGroup}}[(e, e) (e, a) (e, a²) (e, a³); (a, e) (a, a) (a, a²) (a, a³)]
    @show iter = elements(ProductGroup(D6, Z2))
    @show collect(iter) # Tuple{GroupElement{DihedralGroup}, GroupElement{CyclicGroup}}[(e, e) (e, a); (r, e) (r, a); (r², e) (r², a); (s, e) (s, a); (sr, e) (sr, a); (sr², e) (sr², a)]
    @show identity_element(Z2xZ4) # (e, e)
    @show identity_element(ProductGroup(D6, Z2)) # (e, e)
    @show inverse(ae) # (a³, e)
    @show inverse(sa) # (s, a³)
    @show ae * inverse(ae) # (e, e)
    @show sa * inverse(sa) # (e, e)
    @show sa * sa # (s, a) * (s, a) = (s*s, a*a) = (e, a²)
    @show sa * sa == GroupElement((e, aZ2*aZ2), ProductGroup(D6, Z2)) # true
    @show Z2Z2Z2 = ProductGroup(Z2, Z2, Z2) # ℤ₂ × ℤ₂ × ℤ₂
    @show identity_element(Z2Z2Z2) # (e, e, e)
    # @show inverse(GroupElement((a, a, a), Z2Z2Z2)) ERROR: Element a is not in group ℤ₂
    
    @show inverse(GroupElement((aZ2, aZ2, aZ2), Z2Z2Z2)) # (a, a, a)
    aaa = GroupElement((aZ2, aZ2, aZ2), Z2Z2Z2)
    eee = identity_element(Z2Z2Z2)
    @show aaa * aaa # (a, a, a) * (a, a, a) = (a*a, a*a, a*a) = (e, e, e)
    @show aaa * aaa == identity_element(Z2Z2Z2) # true
    @show multiply((aaa, aaa, eee, aaa)) # (a, a, a)
    Z2Z2 = ProductGroup(Z2, Z2)
    @show x = elements(Z2Z2)
    @show ee = x[1] # (e, e)
    @show iter = collect(group_tree(ee, 2))
    @show collect(group_tree(eee, 2))

    @show Z = IntegerGroup() # ℤ
    @show GroupElement(3, Z) # 𝟙³
end

function test_obj()
    D6 = DihedralGroup(3)
    e = identity_element(D6)
    s = GroupElement((1,0), D6)
    r = GroupElement((0,1), D6)
    @show A = Obj(e=>2, s=>3, r=>2, s*r=>15) # 2r ⊕ 15sr ⊕ 3s ⊕ 2e
    @show zero_obj(D6) # 0e
    @show dual_obj(A) # 3s ⊕ 2r² ⊕ 2e ⊕ 15sr
    @show A == dual_obj(dual_obj(A)) # true

    Z = IntegerGroup()
    zero = identity_element(Z)
    one = GroupElement(1, Z)
    @show one * one
    @show inverse(one)
    @show zero * one == one # true
    @show zero * zero == zero # true
    @show B = Obj(zero=>2, one=>3) # 2𝟙⁰ ⊕ 3𝟙¹
    @show dual_obj(B) # 2𝟙⁰ ⊕ 3𝟙⁻¹
    @show B == dual_obj(dual_obj(B)) # true
    @show zero_obj(Z) # 0𝟙⁰

    Z2 = CyclicGroup(2)
    e = identity_element(Z2)
    a = GroupElement(1, Z2)
    @show C = Obj(e=>2, a=>3) # 2e ⊕ 3a
    @show zero_obj(Z2) # 0e
    Z2Z2 = ProductGroup(Z2, Z2)
    ee = identity_element(Z2Z2)
    ea = GroupElement((e, a), Z2Z2)
    ae = GroupElement((a, e), Z2Z2)
    aa = GroupElement((a, a), Z2Z2)
    @show D = Obj(ee=>2, ea=>3, ae=>4, aa=>5) # 3(e, a) ⊕ 2(e, e) ⊕ 5(a, a) ⊕ 4(a, e)
    @show zero_obj(Z2Z2) # 0(e, e)
    @show dual_obj(D) # 3(e, a) ⊕ 4(a, e) ⊕ 5(a, a) ⊕ 2(e, e)

    Z3 = CyclicGroup(3)
    e = identity_element(Z3)
    a = GroupElement(1, Z3)
    b = GroupElement(2, Z3)
    @show E = Obj(e=>2, a=>3, b=>4) # 2e ⊕ 3a ⊕ 4a²
    @show zero_obj(Z3) # 0e
    @show dual_obj(E) # 3a² ⊕ 4a ⊕ 2e
    @show E == dual_obj(dual_obj(E)) # true

    Z2Z = ProductGroup(Z2, Z)
    e = identity_element(Z2Z)
    a = GroupElement(1, Z2)
    one = GroupElement(1, Z)
    @show a1 = GroupElement((a, one), Z2Z) # (a, 𝟙¹)
    @show a1 * a1 # (e, 𝟙²)
    @show inverse(a1) # (a, 𝟙⁻¹)
    @show zero_obj(Z2Z) # (e, 𝟙⁰)
    @show F = Obj(e=>2, a1=>3) # 2(e, 𝟙⁰) ⊕ 3(a, 𝟙¹)
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

    Z2 = CyclicGroup(2)
    e = identity_element(Z2)
    a = GroupElement(1, Z2)
    @show sect = Sector(a, a)
    @show sect[1]
    @show sect[2]

    Z2Z2 = ProductGroup(Z2, Z2)
    ee = identity_element(Z2Z2)
    ea = GroupElement((e, a), Z2Z2)
    ae = GroupElement((a, e), Z2Z2)
    aa = GroupElement((a, a), Z2Z2)
    @show sect = Sector(ee, ea, ae, aa)
    @show sect[1]
    @show sect[2]

    Z3 = CyclicGroup(3)
    e = identity_element(Z3)
    a = GroupElement(1, Z3)
    b = GroupElement(2, Z3)
    @show sect = Sector(a, e, b) # a ⊗ e ⊗ a²

    Z = IntegerGroup()
    Z2Z = ProductGroup(Z2, Z)
    e = identity_element(Z2Z)
    a = GroupElement(1, Z2)
    one = GroupElement(1, Z)
    a1 = GroupElement((a, one), Z2Z)
    @show sect = Sector(e, a1, inverse(a1)) # (e, 𝟙⁰) ⊗ (a, 𝟙¹) ⊗ (a, 𝟙⁻¹)

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
    # @show T # The data is too long to show. The e ⊗ e ⊗ e ⊗ e sector is rand(2,2,2,2) given above
    # @show T[e,s,s,e]  = rand(2,4,3,2)
    # @show T # The data is too long to show. The e ⊗ s ⊗ s ⊗ e sector is rand(2,4,3,2) given above
    @show TT = random_mor(Float64, (A, B, C, D)) # The data is too long to show. Random morphism is generated   
    # @show TT = zero_mor(Float64, (A, B, C, D)) # The data is too long to show. Zero morphism is generated
    # @show Mat = identity_mor(Float64, A) # The data is too long to show. Identity morphism is generated

    
end

function test_Mor_Z2()
    Z2 = CyclicGroup(2)
    e = identity_element(Z2)
    a = GroupElement(1, Z2)
    B = Obj(e=>1, a=>1)
    C = Obj(e=>1, a=>1)
    @show T = random_mor(Float64, (B,C))
    @show zero_mor(Float64, (B,C))
    A = Obj(e=>2)
    @show Mat = identity_mor(Float64, A) # Mor{CyclicGroup, Float64}((2e, 2e), Dict{Sector{CyclicGroup}, Array{Float64}}(e ⊗ e => [1.0 0.0; 0.0 1.0], a ⊗ a => Matrix{Float64}(undef, 0, 0)))
    @show Mat2 = identity_mor(Float64, B) # Mor{CyclicGroup, Float64}((e ⊕ a, e ⊕ a), Dict{Sector{CyclicGroup}, Array{Float64}}(e ⊗ e => [1.0;;], a ⊗ a => [1.0;;]))
end

function test_Mor_D6()
    D6 = DihedralGroup(3)
    e = identity_element(D6)
    s = GroupElement((1,0), D6)
    r = GroupElement((0,1), D6)
    A = Obj(e=>1, s=>1, r=>2, s*r=>1)
    B = Obj(e=>1, s=>1, r*r=>1, s*r=>1)
    @show T = random_mor(Float64, (A,B))
    @show zero_mor(Float64, (A,B))
    @show identity_mor(Float64, A)
end

function test_get_indices()
    D6 = DihedralGroup(3)
    e = identity_element(D6)
    s = GroupElement((1,0), D6)
    r = GroupElement((0,1), D6)
    A = Obj(e=>2, s=>3, r=>2, s*r=>2) 
    B = Obj(e=>2, s=>2, r=>1, s*r=>2)
    C = Obj(e=>2, s=>3, r=>2, s*r=>2)
    D = Obj(e=>2, s=>2, r=>2, s*r*r=>1)
    T = random_mor(Float64, (A, B, C, D))

    @show T[e,e,e,e] # 2×2×2×2 Array{Float64,4}
    @show T[1] # 2e ⊕ 2r ⊕ 3s ⊕ 2sr
    @show T[2] # 2e ⊕ r ⊕ 2s ⊕ 2sr
    @show T[2:3] # (2e ⊕ r ⊕ 2s ⊕ 2sr, 2e ⊕ 2r ⊕ 3s ⊕ 2sr)
    @show T[end] # 2e ⊕ 2r ⊕ 2s ⊕ sr²
    @show lastindex(T) # 4
    @show get_sector_size(T, (e,e,e,e)) # (2, 2, 2, 2)
    @show size = get_sector_size(T, (e,s,s,e)) # (2, 2, 3, 2)
    # @show get_sector_size(T, (e,s,r,e)) ERROR: ArgumentError: The sector (e, s, r, e) is not consistent
    S = Sector(e, s, s, e)
    @show T[S] = rand(size...) # The data is too long to show. Random data is generated and assigned to the sector (e, s, s, e)

end

function test_accend()
    @show is_accend((4,1), 4)
    @show is_accend((4,1,2), 4)
    @show is_accend((4,1,2,3), 4)
    @show is_accend((4,1,2,3,4), 4)
    @show is_accend((1,2), 4)
    @show is_accend((1,2,3), 4)
    @show is_accend((1,2,3,4), 4)
    @show is_accend((1,2,4,1), 4) # false
    @show is_descend((4,1), 4) # false
    @show is_descend((3,2,1), 4) 
    @show is_descend((2,1,4), 4)
    @show is_descend((4,3,2), 4)

    @show is_cyclic((4,1,2,3), 4) # true
    @show is_cyclic((4,1,2,4), 4) # false
    @show is_cyclic((4,1,2,3,4), 4) # false

    @show to_perm((1,2), 4)
    @show to_perm((1,2,3), 4)
    @show to_perm((1,2,3,4), 4)
    @show to_perm((4,1), 4)
    @show to_perm((4,1,2), 4)
    @show to_perm((4,1,2,3), 4)
    @show to_perm((3,4,1), 4)
    # @show to_perm((4,1,2,3,4), 4) # ERROR: ArgumentError: The length of the tuple (4, 1, 2, 3, 4) is greater than 4
    # @show to_perm((1,2,4,1), 4) # ERROR: ArgumentError: The (1, 2, 4, 1) is not in an accending order
end

function test_permute()
    Z3 = CyclicGroup(3)
    e = identity_element(Z3)
    a = GroupElement(1, Z3)
    b = GroupElement(2, Z3)
    S = Sector(a, e, b)
    @show VecG_permutesectors(S, (1,2,3)) # a ⊗ e ⊗ a²
    @show VecG_permutesectors(S, (2,3,1)) # e ⊗ a² ⊗ a
    @show VecG_permutesectors(S, (3,1,2)) # a² ⊗ a ⊗ e

    D6 = DihedralGroup(3)
    e = identity_element(D6)
    s = GroupElement((1,0), D6)
    r = GroupElement((0,1), D6)
    S = Sector(e, s, r, s*r)
    @show VecG_permutesectors(S, (1,2,3,4)) # e ⊗ s ⊗ r ⊗ sr
    @show VecG_permutesectors(S, (2,3,4,1)) # s ⊗ r ⊗ sr ⊗ e
    @show VecG_permutesectors(S, (3,4,1,2)) # r ⊗ sr ⊗ e ⊗ s

    A = Obj(e=>2, s=>3, r=>2, s*r=>1)
    B = Obj(e=>2, s=>3, r=>2, s*r=>1)
    C = Obj(e=>2, s=>3, r=>2, s*r=>1)
    D = Obj(e=>2, s=>3, r=>2, s*r=>1)
    T = random_mor(Float64, (A, B, C, D))
    Tperm = VecG_permutedims(T, (2,3,4,1)) 
    @show Tperm[e,s,r,s*r] == permutedims(T[s*r, e, s, r],(2,3,4,1)) # true

    Z2 = CyclicGroup(2)
    e = identity_element(Z2)
    a = GroupElement(1, Z2)
    B = Obj(e=>2, a=>2)
    C = Obj(e=>2, a=>2)
    @show T = random_mor(Float64, (B,C))
    @show Tperm = VecG_permutedims(T, (2,1))
end

function test_add()
    G = CyclicGroup(3)
    e = GroupElement(0, G)
    a = GroupElement(1, G)
    aa = GroupElement(2, G)
    A = Obj(e=>1, a=>2, aa=>3)
    B = Obj(e=>2, a=>3, aa=>2)
    C = Obj(e=>1, a=>2, aa=>3)
    D = Obj(e=>2, a=>3, aa=>2)
    T1 = random_mor(Float64, (A, B, C, D))
    T2 = random_mor(Float64, (A, B, C, D))
    Tadd = T1 + T2
    @show Tadd[e,e,e,e] == T1[e,e,e,e] + T2[e,e,e,e] # true
    @show Tadd[e,a,a,a] == T1[e,a,a,a] + T2[e,a,a,a] # true
    @show Tadd[a,a,e,a] == T1[a,a,e,a] + T2[a,a,e,a] # true

    Tminus = T1 - T2
    @show Tminus[e,e,e,e] == T1[e,e,e,e] - T2[e,e,e,e] # true
    @show Tminus[e,a,a,a] == T1[e,a,a,a] - T2[e,a,a,a] # true
    @show Tminus[a,a,e,a] == T1[a,a,e,a] - T2[a,a,e,a] # true

    Tdiv = T1 ./ 2
    @show Tdiv[e,e,e,e] == T1[e,e,e,e] ./ 2 # true
    @show Tdiv[e,a,a,a] == T1[e,a,a,a] ./ 2 # true
    @show Tdiv[a,a,e,a] == T1[a,a,e,a] ./ 2 # true

    Tsqrt = sqrt.(T1)
    @show Tsqrt[e,e,e,e] == sqrt.(T1[e,e,e,e]) # true
    
end

function test_dag()
    G = CyclicGroup(2)
    e = GroupElement(0, G)
    a = GroupElement(1, G)
    A = Obj(e=>1, a=>2)
    B = Obj(e=>2, a=>1)
    T = random_mor(ComplexF64, (A, B))
    @show Tdag = VecG_dag(T)
    @show Tdag[e,e] == T[e,e]' # true
    @show [1 + im 2 + 3im]' 
    @show Tdag[a,a] == T[a,a]' # true
    @show T[e,e] = reshape(Array{Complex{Float64}, 2}([1.0+2im 2+3.0im]), 1, 2)
    T[a,a] = reshape(Array{Complex{Float64}, 2}([1.0+2im 2+3.0im]), 2, 1)
    Tdag = VecG_dag(T)
    @show Tdag[e,e]
    @show Tdag[a,a]
    @show Tdag[a,a] == T[a,a]' # true

    G = CyclicGroup(3)
    e = GroupElement(0, G)
    a = GroupElement(1, G)
    aa = GroupElement(2, G)
    A = Obj(e=>1, a=>2, aa=>3)
    B = Obj(e=>2, a=>4, aa=>2)
    @show T = random_mor(ComplexF64, (A, B))
    @show Tdag = VecG_dag(T)

    G = DihedralGroup(3)
    e = identity_element(G)
    s = GroupElement((1,0), G)
    r = GroupElement((0,1), G)
    A = Obj(e=>1, s=>2, r=>3, s*r=>4)
    B = Obj(e=>2, s=>3, r=>2, s*r=>1)
    C = Obj(e=>1, s=>2, r=>3, s*r=>4)
    D = Obj(e=>2, s=>3, r=>2, s*r=>1)
    @show T = random_mor(ComplexF64, (A, B, C, D))
    @show T[e,s,r,s*r]
    @show Tdag = VecG_dag(T)
    @show keys(Tdag.data)
    @show Tdag[e,e,e,e]
    @show Tdag.objects
    @show Tdag[e, e, s*r, s*r] == conj.(permutedims(T[s*r, s*r, e, e], (4,3,2,1))) # true

end

function test_maxabs()
    G = CyclicGroup(2)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       A = Obj(e=>1, a=>1)
       B = Obj(e=>1, a=>1)
         T = Mor(Float64, (A, B))
         T[e,e] = Array{Float64, 2}(reshape([1],1,1)) # Explicit type conversion is needed. Otherwise, the type of the array will be an error
            T[a,a] =Array{Float64, 2}(reshape([2],1,1))
       @show max_abs(T) == 2.0 # true

         G = CyclicGroup(3)
         e = GroupElement(0, G)
            a = GroupElement(1, G)
            aa = GroupElement(2, G)
            A = Obj(e=>1, a=>2, aa=>3)
            B = Obj(e=>1, a=>2, aa=>3)
            T = random_mor(Float64, (A, B))
            @show max_abs(T)
end
# test_group()
# test_obj()
# test_sector()
# test_Mor()
# test_Mor_Z2()
# test_Mor_D6()
# test_get_indices()
# test_accend()
# test_permute()
# test_add()
# test_dag()
test_maxabs()
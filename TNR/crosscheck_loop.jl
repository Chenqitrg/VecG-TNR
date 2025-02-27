using Revise
includet("VecGloop.jl")
using .VecG_TNR
using TensorOperations

function test_dag()
    Z2(i::Int) = GroupElement(i, CyclicGroup(2))
    e = Z2(0)
    a = Z2(1)
    # a2 = Z3(2)
    I = Obj(e=>1)
    J = Obj(e=>1)
    K = Obj(e=>3, a=>3)
    L = Obj(e=>2, a=>1)
    morA = random_mor(ComplexF64, (I, J, K, L))
    dagA = VecG_dag(morA)
    # SS = VecG_tensordot(S, VecG_dag(S), (2,), (2,))
    TT = VecG_tensordot(morA, dagA, 2)
    # @show reshape(morA[e,e,e,e],1, 6) * reshape(dagA[e,e,e,e], 6,1)
    @show TT[e,e,e,e]
    # return morA, dagA
    # @show SS[e,e,e,e]
end

function test_dag0()
    A = rand(ComplexF64, 2,3,4)
    B = permutedims(conj.(A), (3,2,1))
    @show reshape(A,2,12) * reshape(B, 12,2)
    return A, B
end

function test_TT()
    Z3(i::Int) = GroupElement(i, CyclicGroup(3))
    e = Z3(0)
    a = Z3(1)
    I = Obj(e=>1, a=>1, a*a=>2)
    J = Obj(e=>2, a=>1, a*a=>3)
    K = Obj(e=>3, a=>2, a*a=>1)
    L = Obj(e=>4, a=>1, a*a=>4)
    morA = random_mor(ComplexF64, (I, J, K, L))
    morB = random_mor(ComplexF64, (dual_obj(K), dual_obj(L), dual_obj(I), dual_obj(J)))

    loopT = to_T_array(morA, morB)
    loopTT = to_TT_array(loopT)
    loopS = to_S_array(loopT, 33, 1e-15)
    loopSS = to_SS_array(loopS)
    loopTSS = to_TSS_array(loopT, loopS)

    @show to_number(loopTT)
    @show to_number(loopSS)
    @show to_number(loopTSS)
    @show cost = cost_function(loopTT, loopTSS, loopSS)
    # return loopTT, cost

end


function test_to_number()
    Z2(i::Int) = GroupElement(i, CyclicGroup(2))
    e = Z2(0)
    a = Z2(1)
    I = Obj(e=>1, a=>1)
    J = Obj(e=>1, a=>1)
    K = Obj(e=>1, a=>1)
    L = Obj(e=>1, a=>1)
    morA = random_mor(ComplexF64, (I, J, K, L))
    morB = random_mor(ComplexF64, (dual_obj(K), dual_obj(L), dual_obj(I), dual_obj(J)))

    loopT = to_T_array(morA, morB)
    loopTT = to_TT_array(loopT)
    loopS = to_S_array(loopT, 10, 1e-15)
    loopSS = to_SS_array(loopS)
    loopTSS = to_TSS_array(loopT, loopS)

    T = VecG_tensordot(loopT[1], loopT[2], (4,), (1,))
    T = VecG_tensordot(T, loopT[3], (6,), (1,))
    T = VecG_tensordot(T, loopT[4], (8,), (1,))
    T = VecG_permutedims(T, (2,3,4,5,6,7,8,9,10,1))
    T = VecG_partial_trace(T, 1)

    TT = VecG_square(T)
    TTp = to_number(loopTT)
    return TT, TTp
end

function test_to_N()
    group = CyclicGroup(3)
    Z3(i::Int) = GroupElement(i, CyclicGroup(3))
    e = Z3(0)
    a = Z3(1)
    I = Obj(e=>1, a=>1, a*a=>2)
    J = Obj(e=>2, a=>1, a*a=>3)
    K = Obj(e=>3, a=>2, a*a=>1)
    L = Obj(e=>4, a=>1, a*a=>4)
    morA = random_mor(ComplexF64, (I, J, K, L))
    morB = random_mor(ComplexF64, (dual_obj(K), dual_obj(L), dual_obj(I), dual_obj(J)))

    loopT = to_T_array(morA, morB)
    loopTT = to_TT_array(loopT)
    loopS = to_S_array(loopT, 33, 1e-15)
    loopSS = to_SS_array(loopS)
    loopTSS = to_TSS_array(loopT, loopS)

    site = 8
    N = to_N(loopS, loopSS, site)
    number = 0.
    for g_circ in elements(group)
        for insect in group_tree(e, 3)
            outsect = (g_circ * insect[1], insect[2], insect[3]*inverse(g_circ))
            tubesect = TubeSector(insect, g_circ, outsect)
            N_sect = N[tubesect]
            S = loopS[site]
            S_sect = S[insect...]
            S_sect_dag = VecG_dag(S)[(reverse(inverse.(outsect)))...]
            @show tubesect
            @show size(N_sect)
            @show size(S_sect)
            @show size(S_sect_dag)
            @tensor temp_numb[] := N_sect[x,y,z, u,v,w] * S_sect[x,y,z] * S_sect_dag[w,v,u]
            number = number + temp_numb[1]
        end
    end


    return number, to_number(loopSS)
end

function test_to_W()
    group = CyclicGroup(3)
    Z3(i::Int) = GroupElement(i, CyclicGroup(3))
    e = Z3(0)
    a = Z3(1)
    I = Obj(e=>1, a=>1, a*a=>2)
    J = Obj(e=>2, a=>1, a*a=>3)
    K = Obj(e=>3, a=>2, a*a=>1)
    L = Obj(e=>4, a=>1, a*a=>4)
    morA = random_mor(ComplexF64, (I, J, K, L))
    morB = random_mor(ComplexF64, (dual_obj(K), dual_obj(L), dual_obj(I), dual_obj(J)))

    loopT = to_T_array(morA, morB)
    loopTT = to_TT_array(loopT)
    loopS = to_S_array(loopT, 33, 1e-15)
    loopSS = to_SS_array(loopS)
    loopTSS = to_TSS_array(loopT, loopS)

    site =1
    N = to_N(loopS, loopSS, site)
    W = to_W(loopT, loopS, loopTSS, site)
    id = Obj(e=>1)

    W_p = zero_mor(ComplexF64, (Obj(e=>1),), W.out)
    for g_circ in elements(group)
        for insect in group_tree(e, 3)
            outsect = (g_circ * insect[1], insect[2], insect[3]*inverse(g_circ))
            tubesect = TubeSector(insect, g_circ, outsect)
            Wp_sect = TubeSector((e,), g_circ, outsect)
            N_sect = N[tubesect]
            S = loopS[site]
            S_sect = S[insect...]
            @tensor temp_numb[u,v,w] := N_sect[x,y,z, u,v,w] * S_sect[x,y,z]
            W_p[Wp_sect] = W_p[Wp_sect] + reshape(temp_numb, (1, size(temp_numb)...))
        end
    end

    b = (W[(e,),e,(e,e,e)] ≈ W_p[(e,),e,(e,e,e)])
    
    return b
end




Z3(i::Int) = GroupElement(i, CyclicGroup(3))
e = Z3(0)
a = Z3(1)
test_to_W()


# Z3(i::Int) = GroupElement(i, CyclicGroup(3))
#     e = Z3(0)
#     a = Z3(1)
#     a2 = Z3(2)
# loopTT, cost = test_TT()

# A = VecG_permutedims(loopTT[1], (4,1,2,3))
# morA[e,e,e,e]
# dagA[e,e,e,e]

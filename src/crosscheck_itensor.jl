using ITensors
using Revise
includet("main.jl")
using .VecG_TNR


function check_Z2_qr()
    q0 = QN(0,2)
    q1 = QN(1,2)
    i = Index(q0=>1, q1=>1)
    j = Index(q0=>1, q1=>1)
    k = Index(q0=>1, q1=>1)
    l = Index(q0=>1, q1=>1)
    T = randomITensor(i,j,k,l)

    Z2(i::Int) = GroupElement(i-1, CyclicGroup(2))
    e = Z2(1)
    a = Z2(2)
    I = Obj(e=>1, a=>1)
    J = Obj(e=>1, a=>1)
    K = Obj(e=>1, a=>1)
    L = Obj(e=>1, a=>1)
    mor = Mor(Float64, (I, J, K, L))

    for x = 1 : 2, y = 1 : 2, z = 1 : 2
        w = mod((x-1) + (y-1) + (z - 1) ,2) + 1
        mor[Z2(x), Z2(y), Z2(z), Z2(w)] = reshape([T[i=>x, j=>y, k=>z, l=>w]], 1,1,1,1)
    end

    q, r = qr(T, (i, j))
    Q, R = VecG_qr(mor, 2)

    return q, r, Q, R
end

function check_Z3_svd()
    q0 = QN(0,3)
    q1 = QN(1,3)
    q2 = QN(2,3)
    i = Index(q0=>1, q1=>2, q2=>3)
    j = Index(q0=>2, q1=>2, q2=>5)
    k = Index(q0=>1, q1=>3, q2=>3)
    l = Index(q0=>1, q1=>1, q2=>3)
    T = randomITensor(ComplexF64, i,j,k,l)

    arr = array(T, i, j, k, l)

    Z3(i::Int) = GroupElement(i, CyclicGroup(3))
    e = Z3(0)
    a = Z3(1)
    a2 = Z3(2)
    I = Obj(e=>1, a=>2, a2=>3)
    J = Obj(e=>2, a=>2, a2=>5)
    K = Obj(e=>1, a=>3, a2=>3)
    L = Obj(e=>1, a=>1, a2=>3)
    mor = Mor(ComplexF64, (I, J, K, L))

    range1 = Dict(e=>1:1, a=>2:3, a2=>4:6)
    range2 = Dict(e=>1:2, a=>3:4, a2=>5:9)
    range3 = Dict(e=>1:1, a=>2:4, a2=>5:7)
    range4 = Dict(e=>1:1, a=>2:2, a2=>3:5)

    for tup in group_tree(e, 4)
        mor[tup...] = arr[range1[tup[1]], range2[tup[2]], range3[tup[3]], range4[tup[4]]]
    end

    u,s,v = svd(T, (i,j))
    U, S, V = VecG_svd(mor, (1,2), 5)

    return bo
end

function check_Z3_svd(T::Type)
    Z3(i::Int) = GroupElement(i, CyclicGroup(3))
    e = Z3(0)
    a = Z3(1)
    a2 = Z3(2)
    I = Obj(e=>1, a=>2, a2=>3)
    J = Obj(e=>2, a=>2, a2=>5)
    K = Obj(e=>1, a=>3, a2=>3)
    L = Obj(e=>1, a=>1, a2=>3)
    mor = random_mor(T, (I, J, K, L))

    U, S, V = VecG_svd(mor, (1,2), 100)
    T1 = VecG_tensordot(U, S, (3,),(1,))
    morp = VecG_tensordot(T1, V, (3,), (3,))
    return mor, morp
end

function check_Z3_qr()
    q0 = QN(0,3)
    q1 = QN(1,3)
    q2 = QN(2,3)
    i = Index(q0=>1, q1=>2, q2=>3)
    j = Index(q0=>2, q1=>2, q2=>5)
    k = Index(q0=>1, q1=>3, q2=>3)
    l = Index(q0=>1, q1=>1, q2=>3)
    T = randomITensor(i,j,k,l)

    arr = array(T, i, j, k, l)

    Z3(i::Int) = GroupElement(i, CyclicGroup(3))
    e = Z3(0)
    a = Z3(1)
    a2 = Z3(2)
    I = Obj(e=>1, a=>2, a2=>3)
    J = Obj(e=>2, a=>2, a2=>5)
    K = Obj(e=>1, a=>3, a2=>3)
    L = Obj(e=>1, a=>1, a2=>3)
    mor = Mor(Float64, (I, J, K, L))

    range1 = Dict(e=>1:1, a=>2:3, a2=>4:6)
    range2 = Dict(e=>1:2, a=>3:4, a2=>5:9)
    range3 = Dict(e=>1:1, a=>2:4, a2=>5:7)
    range4 = Dict(e=>1:1, a=>2:2, a2=>3:5)

    for tup in group_tree(e, 4)
        mor[tup...] = arr[range1[tup[1]], range2[tup[2]], range3[tup[3]], range4[tup[4]]]
    end

    q,r = qr(T, (i,j))
    Q,R = VecG_qr(mor, 2)
    return q,r,Q,R # Result: not 100% match, just some numbers match.
end

function check_Z3_contract()
    q0 = QN(0,3)
    q1 = QN(1,3)
    q2 = QN(2,3)
    i = Index(q0=>1, q1=>2, q2=>3)
    j = Index(q0=>2, q1=>2, q2=>5)
    k = Index(q0=>1, q1=>3, q2=>3)
    l = Index(q0=>1, q1=>1, q2=>3)
    T = randomITensor(i,j,k,l)

    arr = array(T, i, j, k, l)

    Z3(i::Int) = GroupElement(i, CyclicGroup(3))
    e = Z3(0)
    a = Z3(1)
    a2 = Z3(2)
    I = Obj(e=>1, a=>2, a2=>3)
    J = Obj(e=>2, a=>2, a2=>5)
    K = Obj(e=>1, a=>3, a2=>3)
    L = Obj(e=>1, a=>1, a2=>3)
    mor = Mor(Float64, (I, J, K, L))

    range1 = Dict(e=>1:1, a=>2:3, a2=>4:6)
    range2 = Dict(e=>1:2, a=>3:4, a2=>5:9)
    range3 = Dict(e=>1:1, a=>2:4, a2=>5:7)
    range4 = Dict(e=>1:1, a=>2:2, a2=>3:5)

    for tup in group_tree(e, 4)
        mor[tup...] = arr[range1[tup[1]], range2[tup[2]], range3[tup[3]], range4[tup[4]]]
    end

    Q,R = VecG_qr(mor, 2)
    mor_p = VecG_tensordot(Q, R, 1)
    return mor, mor_p
end

function check_Z3_contract2()
    q0 = QN(0,3)
    q1 = QN(1,3)
    q2 = QN(2,3)
    i = Index(q0=>1, q1=>2, q2=>3)
    j = Index(q0=>2, q1=>2, q2=>2)
    k = Index(q0=>3, q1=>3, q2=>3)
    l = Index(q0=>4, q1=>1, q2=>3)
    T = randomITensor(i,j,k,l)

    arr = array(T, i, j, k, l)

    Z3(i::Int) = GroupElement(i, CyclicGroup(3))
    e = Z3(0)
    a = Z3(1)
    a2 = Z3(2)
    I = Obj(e=>1, a=>2, a2=>3)
    J = Obj(e=>2, a=>2, a2=>2)
    K = Obj(e=>3, a=>3, a2=>3)
    L = Obj(e=>4, a=>1, a2=>3)
    mor = Mor(Float64, (I, J, K, L))

    range1 = Dict(e=>1:1, a=>2:3, a2=>4:6)
    range2 = Dict(e=>1:2, a=>3:4, a2=>5:6)
    range3 = Dict(e=>1:3, a=>4:6, a2=>7:9)
    range4 = Dict(e=>1:4, a=>5:5, a2=>6:8)

    for tup in group_tree(e, 4)
        mor[tup...] = arr[range1[tup[1]], range2[tup[2]], range3[tup[3]], range4[tup[4]]]
    end

    TT = T * dag(δ(k,dag(k'))) * dag(δ(l,dag(l'))) * dag(T')

    mm = VecG_tensordot(mor, VecG_dag(mor), (3,4), (2,1))
    return TT, mm
end

function check_partial_trace()
    Z3(i::Int) = GroupElement(i, CyclicGroup(3))
    e = Z3(0)
    a = Z3(1)
    a2 = Z3(2)
    I = Obj(e=>1, a=>2, a2=>3)
    J = Obj(e=>2, a=>2, a2=>5)
    K = Obj(e=>1, a=>2, a2=>3)
    L = Obj(e=>2, a=>4, a2=>5)
    M = Obj(e=>2, a=>5, a2=>4)
    N = Obj(e=>1, a=>3, a2=>2)
    mor = random_mor(Float64, (I, J, K, L, M, N))

    newmor = VecG_partial_trace(mor, 2)
end

function check_partial_trace_full()
    Z3(i::Int) = GroupElement(i, CyclicGroup(3))
    e = Z3(0)
    a = Z3(1)
    a2 = Z3(2)
    K = Obj(e=>1, a=>2, a2=>3)
    L = Obj(e=>2, a=>4, a2=>5)
    M = Obj(e=>2, a=>5, a2=>4)
    N = Obj(e=>1, a=>3, a2=>2)
    mor = random_mor(Float64, (K, L, M, N))

    newmor = VecG_partial_trace(mor, 2)
    return newmor
end

function check_add()
    Z3(i::Int) = GroupElement(i, CyclicGroup(3))
    e = Z3(0)
    a = Z3(1)
    a2 = Z3(2)
    I = Obj(e=>1, a=>2, a2=>3)
    J = Obj(e=>2, a=>2, a2=>5)
    K = Obj(e=>1, a=>2, a2=>3)
    L = Obj(e=>2, a=>4, a2=>5)
    M = Obj(e=>2, a=>5, a2=>4)
    N = Obj(e=>1, a=>3, a2=>2)
    mora = random_mor(Float64, (I, J, K, L, M, N))
    morb = random_mor(Float64, (I, J, K, L, M, N))
    return mora, morb, mora + morb, mora - morb
end

function check_square()
    Z2(i::Int) = GroupElement(i, CyclicGroup(1))
    e = Z2(0)
    a = Z2(1)
    I = Obj(e=>1)
    J = Obj(e=>1)
    mora = random_mor(ComplexF64, (I, J))
    @show mora
    @show VecG_dag(mora)
    return VecG_square(mora)
end

function test_tube()
    Z3 = CyclicGroup(3)
    e = GroupElement(0, Z3)
    a = GroupElement(1, Z3)
    aa = GroupElement(2, Z3)
    mor = TubeSector((a,aa,e), a, (aa,a))
    return mor
end

function check_tube_random()
    Z2(i::Int) = GroupElement(i, CyclicGroup(2))
    e = Z2(0)
    a = Z2(1)
    I = Obj(e=>1, a=>1)
    J = Obj(e=>1, a=>1)
    mor = random_mor(Float64, (I,), (J,))
    return mor
end

function check_tube_zero()
    Z2(i::Int) = GroupElement(i, CyclicGroup(2))
    e = Z2(0)
    a = Z2(1)
    I = Obj(e=>1, a=>1)
    J = Obj(e=>1, a=>1)
    mor = zero_mor(Float64, (I,), (J,))
    return mor
end

function check_out_prod()
    S3 = DihedralGroup(3)
    e = GroupElement((0,0), S3)
    s = GroupElement((1,0), S3)
    r = GroupElement((0,1), S3)

    A = Obj(e=>1, s=>2, r=>3, s*r=>2)
    B = Obj(e=>2, s=>2, r=>1, s*(r*r)=>1)
    T1 = random_mor(Float64, (A, A, B))
    T2 = random_mor(Float64, (A, B, B, A))

    return VecG_outprod(T1, T2)
end
# mor = check_partial_trace()
# check_Z3_svd(Float64)
# check_Z3_contract()
# a, b,A, B = check_add()
# check_Z2_svd()
# T1 = VecG_tensordot(U, S, (3,),(1,))
# morp = VecG_tensordot(T1, conj.(V), (3,), (3,))

# check_Z3_contract()
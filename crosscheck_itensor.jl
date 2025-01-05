using ITensors
include("groups.jl")
include("block_matrix_calculation.jl")
include("VecGtensor.jl")
include("display.jl")
include("VecGoperations.jl")


function check_Z2()
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

    f1, k1 = factorize(T, (i, j); which_decomp = "svd", ortho = "none")
    P1, Q1 = VecG_factorize(mor, 2, 4, "svd")

    f2, k2 = factorize(T, (j, k); which_decomp = "svd", ortho = "none")
    perm_mor = VecG_permutedims(mor, (2,3,4,1))
    P2, Q2 = VecG_factorize(perm_mor, 2, 4, "svd")

    return f1, k1, P1, Q1, f2, k2, P2, Q2
end

function check_Z3()
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

    f1, k1 = factorize(T, (i, j); which_decomp = "svd", ortho = "none")
    ind_fact = commonind(f1, k1)
    P1, Q1 = VecG_factorize(mor, 2, 13, "svd")
    return f1, k1, i, j, ind_fact, P1, Q1, i, j
end

Z3(i::Int) = GroupElement(i, CyclicGroup(3))
e = Z3(0)
a = Z3(1)
a2 = Z3(2)
f1, k1, i, j, ind_fact, P1, Q1 = check_Z3()

f1_array = array(f1, i, j,ind_fact)

f1_array[1,1:2,1]
P1[e,e,e]
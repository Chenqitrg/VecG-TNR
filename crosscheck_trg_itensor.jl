include("package/groups.jl")
include("package/VecGtensor.jl")
include("package/block_matrix_calculation.jl")
include("package/VecGdecomposition.jl")
include("package/VecGcontraction.jl")
include("renormalization.jl")
include("package/display.jl")
using ITensors
function init_crosscheck_to_proj(type::Type)
    q0 = QN(0, 3)
    q1 = QN(1, 3)
    q2 = QN(2, 3)
    i = Index(q0 => 1, q1 => 2, q2 => 3;tags = "L,left")
    j = Index(q0 => 2, q1 => 2, q2 => 5; tags = "mid")
    k = Index(q0 => 1, q1 => 3, q2 => 3; tags = "R,right")
    L = randomITensor(type, i, j)
    R = randomITensor(type, k, dag(j))

    L_arr = array(L, i, j)
    R_arr = array(R, k, dag(j))

    Z3(i::Int) = GroupElement(i, CyclicGroup(3))
    e = Z3(0)
    a = Z3(1)
    a2 = Z3(2)
    I = Obj(e => 1, a => 2, a2 => 3)
    J = Obj(e => 2, a => 2, a2 => 5)
    K = Obj(e => 1, a => 3, a2 => 3)
    L_G = Mor(type, (I, J))
    R_G = Mor(type, (K, dual_obj(J)))

    range1 = Dict(e => 1:1, a => 2:3, a2 => 4:6)
    range2 = Dict(e => 1:2, a => 3:4, a2 => 5:9)
    range3 = Dict(e => 1:1, a => 2:4, a2 => 5:7)
    dualrange2 = Dict(e => 1:2, a2 => 3:4, a => 5:9)

    for tup in group_tree(e, 2)
        L_G[tup...] = L_arr[range1[tup[1]], range2[tup[2]]]
        R_G[tup...] = R_arr[range3[tup[1]], dualrange2[tup[2]]]
    end

    # return L, R,L_G, R_G
    return L, R, L_G, R_G
end

function try_dag()
    D4 = DihedralGroup(4)  # 假设这是一个群对象
    e = GroupElement((0, 0), D4)
    s = GroupElement((1, 0), D4)
    r = GroupElement((0, 1), D4)
    r2 = GroupElement((0, 2), D4)
    sr2 = s * r2

    A = Obj(e => 2, s * r => 3, s => 1, r2 => 2, sr2 => 1, r => 5)
    B = Obj(e => 2, s * r => 3, s => 1, r2 => 2, sr2 => 1, r => 5)
    C = Obj(e => 2, s * r => 3, s => 1, r2 => 2, sr2 => 1, r => 5)
    D = Obj(e => 2, s * r => 3, s => 1, r2 => 2, sr2 => 1, r => 5)
    T = random_mor(Float64, (A, B, C, D))

    return VecG_dag(T)
end

function try_single_qr()
    # Mock data
    D4 = DihedralGroup(4)  # 假设这是一个群对象
    e = GroupElement((0, 0), D4)
    s = GroupElement((1, 0), D4)
    r = GroupElement((0, 1), D4)
    r2 = GroupElement((0, 2), D4)
    sr2 = s * r2

    # Create Objects
    A = Obj(e => 2, s * r => 3, s => 1, r2 => 2, sr2 => 1, r => 5)
    B = Obj(e => 2, s * r => 3, s => 1, r2 => 2, sr2 => 1, r => 5)
    C = Obj(e => 2, s * r => 3, s => 1, r2 => 2, sr2 => 1, r => 5)
    # Create BlockTensor
    T = random_mor(Float64, (A, A, A, A))
    L = random_mor(Float64, (A, dual_obj(A)))
    R = random_mor(Float64, (dual_obj(A), A))
    # @show VecG_tensordot(L, R, (2,), (2,))
    return to_projector(L, R, 10^-7)
end

# Given the tag, find the corresponded index 
get_taggedindex(T::ITensor, t::String) = getfirst(x -> hastags(x, t), inds(T))
get_taggedindices(T::ITensor, ts::Tuple{Vararg{String}}) = map(x -> get_taggedindex(T, x), ts)
get_taggedindices(T::ITensor) = get_taggedindices(T, ("left", "up", "right", "down"))

function to_projector(L::ITensor, R::ITensor, epsilon::Float64)
    Ll = get_taggedindex(L, "L, left")
    LR = L * R
    U, S, V = svd(LR, Ll; cutoff=epsilon)
    sqrt_S = sqrt.(S)
    inv_sqrt_S = 1 ./ sqrt_S

    PR = R * dag(V) * dag(inv_sqrt_S)
    PL = dag(inv_sqrt_S) * dag(U) * L

    return PR, PL
end

function try_entanglement_filtering()
    D4 = DihedralGroup(4)  # 假设这是一个群对象
    e = GroupElement((0, 0), D4)
    s = GroupElement((1, 0), D4)
    r = GroupElement((0, 1), D4)
    r2 = GroupElement((0, 2), D4)
    sr2 = s * r2

    A = Obj(e => 2, r => 5, r2 => 2, s => 1, s * r => 3, sr2 => 1)
    B = Obj(e => 2, r => 5, r2 => 2, s => 1, s * r => 4, sr2 => 1)
    C = Obj(e => 2, r => 5, r2 => 3, s => 1, s * r => 3, sr2 => 1)
    D = Obj(e => 2, r => 5, r2 => 2, s => 1, s * r => 3, sr2 => 2)

    morA = random_mor(Float64, (A, C, D, B))
    morB = random_mor(Float64, (dual_obj(D), dual_obj(B), dual_obj(A), dual_obj(C)))

    entanglement_filtering(morA, morB, 10, 1e-7)
end

function try_coase_graining()
    D4 = DihedralGroup(4)  # 假设这是一个群对象
    e = GroupElement((0, 0), D4)
    s = GroupElement((1, 0), D4)
    r = GroupElement((0, 1), D4)
    r2 = GroupElement((0, 2), D4)
    sr2 = s * r2

    A = Obj(e => 2, r => 4, r2 => 2, s => 1, s * r => 3, sr2 => 10)
    B = Obj(e => 2, r => 5, r2 => 2, s => 6, s * r => 4, sr2 => 1)
    C = Obj(e => 2, r => 2, r2 => 3, s => 1, s * r => 3, sr2 => 1)
    D = Obj(e => 2, r => 7, r2 => 2, s => 1, s * r => 3, sr2 => 2)

    morA = random_mor(Float64, (A, C, D, B))
    morB = random_mor(Float64, (dual_obj(D), dual_obj(B), dual_obj(A), dual_obj(C)))

    newmor = coarse_graining(morA, morB, 20)

    return newmor
end

D4 = DihedralGroup(4)  # 假设这是一个群对象
    e = GroupElement((0, 0), D4)
    s = GroupElement((1, 0), D4)
    r = GroupElement((0, 1), D4)
    r2 = GroupElement((0, 2), D4)
    sr2 = s * r2

    try_coase_graining()
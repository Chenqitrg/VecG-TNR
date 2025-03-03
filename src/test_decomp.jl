using Revise
includet("main.jl")
using .VecG_TNR
using ITensors
using LinearAlgebra

function test_sectoutin()
    D6 = DihedralGroup(3)
    e = identity_element(D6)
    s = GroupElement((1,0), D6)
    r = GroupElement((0,1), D6)
    A = Obj(e=>2, s=>3, r=>2, s*r=>1) 
    B = Obj(e=>2, s*r*r=>2, r=>3, s*r=>2)
    C = Obj(e=>2, s=>3, r*r=>2, s*r=>2)
    D = Obj(e=>2, s=>4, r=>1, s*r*r=>1)
    TT = random_mor(Float64, (A, B, C, D)) # The data is too long to show. Random morphism is generated   
    # @show out_sector_tuple, in_sector_tuple = to_sector_outin(keys(TT.data), 2, s)
    # @show to_sector_outin(keys(TT.data), 2, r)
    # @show to_sector_outin(keys(TT.data), 3, s)
    # @show typeof(out_sector_tuple)
    @show sec_split = to_sector_outin(TT, 2)
    Mat = extract_blocks_to_matrix(TT, 2, sec_split)
    @show Mat[s][1,1] == reshape(TT[s,e,r*r,s*r*r], 6,2) # true
    @show Mat[s][1,2] == reshape(TT[s,e,s,e], 6,6) # true
    @show Mat[s][2,2] == reshape(TT[r,s*r,s,e], 4,6) # true
    @show Mat[e][1,1] == reshape(TT[e,e,e,e], 4,4) # true
end
function test_sectoutin_Z2()
    Z2 = CyclicGroup(2)
    e = GroupElement(0, Z2)
    a = GroupElement(1, Z2)
    A = Obj(e=>2, a=>2)
    B = Obj(e=>2, a=>2)
    T = Mor(ComplexF64, (A, B))
    T[e,e] =  Array{ComplexF64}([1+im 2+im; 3-im 4+im])
    T[a,a] =  Array{ComplexF64}([1+im 8+im; 4-im 4+im])
    @show sec_split = to_sector_outin(T, 1) # Dict{Any, Any}(a => (((a,),), ((a,),)), e => (((e,),), ((e,),)))
    @show Mat = extract_blocks_to_matrix(T, 2, sec_split) # Dict{Any, Any}(a => AbstractMatrix[ComplexF64[1.0 + 1.0im; 4.0 - 1.0im; 8.0 + 1.0im; 4.0 + 1.0im;;];;], e => AbstractMatrix[ComplexF64[1.0 + 1.0im; 3.0 - 1.0im; 2.0 + 1.0im; 4.0 + 1.0im;;];;])
end

function test_merge()
    D6 = DihedralGroup(3)
    e = identity_element(D6)
    s = GroupElement((1,0), D6)
    r = GroupElement((0,1), D6)
    A = Obj(e=>2, s=>3, r=>2, s*r=>1) 
    B = Obj(e=>2, s*r*r=>2, r=>3, s*r=>2)
    C = Obj(e=>2, s=>3, r*r=>2, s*r=>2)
    D = Obj(e=>2, s=>4, r=>1, s*r*r=>1)
    TT = random_mor(Float64, (A, B, C, D))
    sec_split = to_sector_outin(TT, 2)
    Mat = extract_blocks_to_matrix(TT, 2, sec_split)
    svd_dic = svd_calculation(Mat)
    for g in keys(svd_dic)
        println("old")
        display(svd_dic[g][2])
    end
    block_cutoff!(svd_dic, 6)
    for g in keys(svd_dic)
        println("new")
        display(svd_dic[g][2])
    end
    # U, S, V = merge_matrix_to_block(TT, 2, svd_dic, sec_split)
    # @show U.objects

end


function check_Z2_svd() # Also a check for svd_dic and merge_matrix_to_block
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

    u, s, v = svd(T, (i, j))

    sec_split = to_sector_outin(mor, 2)
    Mat = extract_blocks_to_matrix(mor, 2, sec_split)
    svd_dic = svd_calculation(Mat)
    U, S, V = merge_matrix_to_block(mor, 2, svd_dic, sec_split)

    @show u
    @show U[e,e,e]
    @show S
    @show s # Two results agrees
end

function check_Z3_svd()
    q0 = QN(0,3)
    q1 = QN(1,3)
    q2 = QN(2,3)
    i = Index(q0=>1, q1=>1, q2=>1)
    j = Index(q0=>1, q1=>1, q2=>1)
    k = Index(q0=>1, q1=>1, q2=>1)
    l = Index(q0=>1, q1=>1, q2=>1)
    T = ITensor(ComplexF64, i,j,k,l)

    # 辅助函数：根据所取块号返回约定的电荷值
# 注意：这里约定
#   块 1 -> 0
#   块 2 -> +1
#   块 3 -> -1
function get_charge(block::Int)
    if block == 1
        return 0
    elseif block == 2
        return 1
    elseif block == 3
        return -1
    else
        error("无效的块号")
    end
end

Matee = zeros(ComplexF64, 3,3)
# 遍历所有组合，只对“中性”（charge neutral）的组合赋非零值
count = 1
for a in 1:3, b in 1:3, c in 1:3, d in 1:3
    # 这里每个数字 a,b,c,d 分别表示在对应 Index 中所取的块号
    if mod(get_charge(a) + get_charge(b) + get_charge(c) + get_charge(d),3) == 0
        # 测试值：可以自定，这里用计数器生成一个复数
        T[i => a, j => b, k => c, l => d] = count + im*count
        count += 1
    end
end
    
# @show T

arr = Array(T, i,j,k,l)

# 第一行
Matee[1,1] = arr[1,1,1,1]  # 行: (1,1), 列: (1,1)
Matee[1,2] = arr[1,1,2,3]  # 行: (1,1), 列: (2,3)
Matee[1,3] = arr[1,1,3,2]  # 行: (1,1), 列: (3,2)

# 第二行
Matee[2,1] = arr[2,3,1,1]  # 行: (2,3), 列: (1,1)
Matee[2,2] = arr[2,3,2,3]  # 行: (2,3), 列: (2,3)
Matee[2,3] = arr[2,3,3,2]  # 行: (2,3), 列: (3,2)

# 第三行
Matee[3,1] = arr[3,2,1,1]  # 行: (3,2), 列: (1,1)
Matee[3,2] = arr[3,2,2,3]  # 行: (3,2), 列: (2,3)
Matee[3,3] = arr[3,2,3,2]  # 行: (3,2), 列: (3,2)

# @show Matee

    Z3(i::Int) = GroupElement(i, CyclicGroup(3))
    e = Z3(0)
    a = Z3(1)
    a2 = Z3(2)
    I = Obj(e=>1, a=>1, a2=>1)
    J = Obj(e=>1, a=>1, a2=>1)
    K = Obj(e=>1, a=>1, a2=>1)
    L = Obj(e=>1, a=>1, a2=>1)
    mor = Mor(ComplexF64, (I, J, K, L))

    range1 = Dict(e=>1:1, a=>2:2, a2=>3:3)
    range2 = Dict(e=>1:1, a=>2:2, a2=>3:3)
    range3 = Dict(e=>1:1, a=>2:2, a2=>3:3)
    range4 = Dict(e=>1:1, a=>2:2, a2=>3:3)

    for tup in group_tree(e, 4)
        mor[tup...] = arr[range1[tup[1]], range2[tup[2]], range3[tup[3]], range4[tup[4]]]
    end

    # @show mor
    # u,s,v = svd(T, i,j; maxdim = 50)
    # u_arr = Array(u, (i, j, commonind(u,s)))
    # s_vec = diag(Array(s, commonind(u,s), commonind(s,v)))
    sec_split = to_sector_outin(mor, 2)
    Mat = extract_blocks_to_matrix(mor, 2, sec_split)
    svd_dic = svd_calculation(Mat)
    U, S, V = VecG_svd(mor, (1,2), 50)
    
    @show sec_split[e]
    Matee_concat = concatenate_matrices_with_metadata(Mat[e])[1]
    # display(Matee_concat)
    # display(Matee)
    ust, sst, vst = svd(Matee)
    ust_con, sst_con, vst_con = svd(Matee_concat)
    # display(ust[1,1:3])
    display(ust_con)
    display(ust)
    display(U[e,e,e])
    
    # display(sst)


    # u, s, v = svd(T, (i, j))

    
    # 
    # U, S, V = merge_matrix_to_block(mor, 2, svd_dic, sec_split)

    # display(Array(u,i,j,commonind(u,s))[1,1,1:3])
    # display(svd_dic[e][1][3])
    # @show sec_split[e][1]
    # @show U[e,e,e]
    # # @show S
    # # @show s # Two results agrees
end

function test_svd_Z2_constructive()
    Z2(i::Int) = GroupElement(i-1, CyclicGroup(2))
    e = Z2(1)
    a = Z2(2)
    I = Obj(e=>3, a=>2)
    J = Obj(e=>4, a=>5)
    T = Mor(ComplexF64, (I,J))
    T[e,e] = Array{ComplexF64, 2}([1 2 3 4; 5 6-im 7 8; 9 10+im 11 12])
    T[a,a] = Array{ComplexF64, 2}([1 2 3+im 4 5;  6+im 7 8 9 10])
    sec_split = to_sector_outin(T,1)
    Mat = extract_blocks_to_matrix(T, 1, sec_split)
    svd_dic = svd_calculation(Mat)
    U, S, V = merge_matrix_to_block(T, 1, svd_dic, sec_split)
    uee, see, vee = svd(T[e,e])
    @show U[e,e] == uee # true
    @show see == diag(S[e,e]) # true
    @show V[e,e] == conj.(vee) # true
end

function test_cutoff()
    Z2(i::Int) = GroupElement(i-1, CyclicGroup(2))
    e = Z2(1)
    a = Z2(2)
    I = Obj(e=>3, a=>2)
    J = Obj(e=>4, a=>5)
    T = Mor(ComplexF64, (I,J))
    T[e,e] = Array{ComplexF64, 2}([1 2 3 4; 5 6-im 7 8; 9 10+im 11 12])
    T[a,a] = Array{ComplexF64, 2}([1 2 3+im 4 5;  6+im 7 8 9 10])
    sec_split = to_sector_outin(T,1)
    Mat = extract_blocks_to_matrix(T, 1, sec_split)
    svd_dic = svd_calculation(Mat)
    Ue, Se, Ve = svd_dic[e]
    display(Ue[1])
    block_cutoff!(svd_dic, 4)
    Ue, Se, Ve = svd_dic[e]
    display(Ue[1])
    # U, S, V = merge_matrix_to_block(T, 1, svd_dic, sec_split)
end

# test_sectoutin()
# test_sectoutin_Z2()
# test_merge()

# check_Z2_svd()
check_Z3_svd()

# test_svd_Z2_constructive()
# test_cutoff()
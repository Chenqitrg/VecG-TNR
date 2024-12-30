include("groups.jl")
include("block_matrix_calculation.jl")
include("VecGtensor.jl")

group = DihedralGroup(4)
legs_to_split = (2,3)

for row_sector in Iterators.product(fill(elements(group), length(legs_to_split))...)
    println(row_sector)
end

function extract_blocks_to_matrix(bt::BlockTensor{Float64, G}, legs_to_split::Vector{Int}, g_bridge::GroupElement{G}) where {G <: Group}
    n = length(bt.objects)

    # Check if the elements of `legs` are within the range 1 to n
    if any(x -> x < 1 || x > n, legs_to_split)
        error("Elements of legs must be in the range 1 to $n")
    end
    
    # Check if `legs` is monotonically increasing and continuous (considering periodic boundary condition)
    for i in 1:length(legs_to_split)-1
        if (legs_to_split[i+1] != legs_to_split[i] + 1) && !(legs_to_split[i] == n && legs_to_split[i+1] == 1)
            error("legs must be monotonically increasing and continuous with periodic boundary condition")
        end
    end

    # Get the rest of the numbers in 1 to n excluding `legs`
    rest = (legs[end]:(legs[end]+n-length(legs)-1))
    rest = map(x -> (x % n) + 1, rest)

    group = bt.group

    for line in Iterators.product(fill(elements(group), length(legs_to_split))...)
        println(line)
    end

    # 假设变量 `vec` 是你要验证的整数向量，变量 `n` 是取值范围的上限
    if isempty(legs_to_split) || !(all(x -> x isa Int && x >= 1 && x <= n, legs_to_split)) || !(issorted(legs_to_split) && all(diff(legs_to_split) .== 1))
        error("The vector is not a sequence of consecutive positive integers within the range 1 to $n.")
    end

    A_blocks = []  # 存储不同 sector 的矩阵块
    A_block_keys = []  # 存储对应的 keys

    # 遍历 BlockTensor 的每个块
    for (key, block) in bt.data
        # 验证是否满足 sector 组合的规则 gigjgk...g_bridge = e
        bridge_compatible = true
        for leg in legs_to_split
            if group.compose(key[leg], g_bridge) != group.identity()
                bridge_compatible = false
                break
            end
        end

        # 如果符合规则，将块数据加入大矩阵 A
        if bridge_compatible
            push!(A_blocks, block)
            push!(A_block_keys, key)
        end
    end

    # 如果没有符合条件的块，返回空
    if isempty(A_blocks)
        return nothing, nothing
    end

    # 合并块为一个大矩阵 A
    A = cat(A_blocks...; dims=1)
    return A, A_block_keys
end


function extract_blocks_to_matrix(bt::BlockTensor{Float64, G}, legs_to_split::Vector{Int}, g_bridge::G) where {G <: GroupElement}
    group = bt.group
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
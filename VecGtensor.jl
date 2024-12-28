include("groups.jl")

# 定义 BlockTensor 结构体
struct BlockTensor{Float64, G <: GroupElement}
    group::Type{G}  # 表示 GroupElement 的具体类型
    dimensions::Vector{Dict{G, Int}}  # 映射 group elements 到维度的向量
    data::Dict{Vector{G}, Array{Float64}}  # 主数据，键是 group elements 的向量，值是数组
end

# 定义构造函数
function BlockTensor(::Type{G}, dimensions::Vector{Dict{G, Int}}) where G <: GroupElement
    data = Dict{Vector{G}, Array{Float64}}()  # 初始化空数据字典
    return BlockTensor{Float64, G}(G, dimensions, data)
end

# Function to set a specific block of the tensor
function set_block!(bt::BlockTensor{Float64, G}, key::Vector{G}, value::Array{Float64}) where {G <: GroupElement}
    # Ensure the key matches the group and dimensions are consistent
    if length(key) != length(bt.dimensions)
        throw(ArgumentError("The key length does not match the number of tensor legs."))
    end
    for (i, g) in enumerate(key)
        if !(g in keys(bt.dimensions[i]))
            throw(ArgumentError("Group element $g is not valid for leg $i."))
        end
        expected_dim = bt.dimensions[i][g]
        if size(value, i) != expected_dim
            throw(ArgumentError("The size of the array along dimension $i does not match the expected dimension $expected_dim for group element $g."))
        end
    end
    bt.data[key] = value
end

# Function to get a block of the tensor
function get_block(bt::BlockTensor{T, G}, key::Vector{G}) where {T, G <: GroupElement}
    return get(bt.data, key, nothing)
end

# Custom display for BlockTensor
function Base.show(io::IO, bt::BlockTensor)
    println(io, "BlockTensor with group: ", bt.group)
    println(io, "Dimensions: ")
    for (i, dim_dict) in enumerate(bt.dimensions)
        println(io, "  Leg $i:")
        for (k, v) in dim_dict
            println(io, "    Group element $k => Dimension $v")
        end
    end
    println(io, "Blocks: ")
    for (key, value) in bt.data
        println(io, "  Key: ", key, " => Block: ", value)
    end
end



# # Example usage:
# # Define a toy group (e.g., integers modulo 3)
# struct ModGroup <: GroupElement
#     n::Int
# end

# Base.:+(g1::ModGroup, g2::ModGroup) = ModGroup((g1.n + g2.n) % 3)
# Base.:-(g1::ModGroup, g2::ModGroup) = ModGroup((g1.n - g2.n) % 3)
# Base.:(==)(g1::ModGroup, g2::ModGroup) = g1.n == g2.n
# Base.hash(g::ModGroup, h::UInt) = hash(g.n, h)
# Base.show(io::IO, g::ModGroup) = print(io, "ModGroup(", g.n, ")")

# # Define dimensions for a 2-leg tensor, where each leg has sectors of specific dimensions
# dim1 = Dict(ModGroup(0) => 2, ModGroup(1) => 3, ModGroup(2) => 1)
# dim2 = Dict(ModGroup(0) => 1, ModGroup(1) => 2)
# dim3 = Dict(ModGroup(0) => 2, ModGroup(1) => 2, ModGroup(2) => 2)

# # Create the BlockTensor
# bt = BlockTensor(ModGroup, [dim1, dim2, dim3])

# set_block!(bt, [ModGroup(0), ModGroup(0), ModGroup(2)], rand(2,1,2))

# # # Set a block
# # set_block!(bt, [ModGroup(0), ModGroup(1)], rand(2, 2))
# # set_block!(bt, [ModGroup(1), ModGroup(0)], rand(3, 1))

# # Get a block
# block = get_block(bt, [ModGroup(0), ModGroup(0), ModGroup(2)])
# println("Retrieved block: ", block)

# # Display the BlockTensor
# # @show bt


include("groups.jl")

# Define the object on the bond
struct Object{G <: Group}
    obj::Dict{GroupElement{G}, Int}
    group::G
    function Object(X::Dict{GroupElement{G}, Int}, group::G) where {G<:Group}
        el = elements(group)
        for g in el
            if !haskey(X, g)
                X[g] = 0
            end
        end
        return new{G}(X, group)
    end
end

# Custom the display form for an object to be direct sum
function Base.show(io::IO, x::Object{G}) where {G<:Group}
    group = x.group
    obj = x.obj
    el = elements(group)
    j = 1
    for (j,i) in enumerate(el)
        if obj[i] != 0
            if j != 1
                print(io, " ⊕ ")
            end
            if obj[i] == 1
                print(io, i)
            else
                print(io, obj[i], i)
            end
        end
    end
end

# Custom the display form of a sector
function Base.show(io::IO, v::Vector{GroupElement{G}}) where {G<:Group}
    for (i, g) in enumerate(v)
        if i != 1
            print(io, " ⊗ ")
        end
        print(io, g)
    end
end

# 定义 BlockTensor 结构体
struct BlockTensor{Float64, G <: Group}
    group::G  # 表示 GroupElement 的具体类型
    objects::Vector{Object{G}}  # 映射 group elements 到维度的向量
    data::Dict{Vector{GroupElement{G}}, Array{Float64}}  # 主数据，键是 group elements 的向量，值是数组
end


# 定义构造函数
function BlockTensor(group::G, objects::Vector{Object{G}}) where {G<:Group}
    data = Dict{Vector{GroupElement{G}}, Array{Float64}}()  # 初始化空数据字典
    return BlockTensor{Float64, G}(group, objects, data)
end

# Custom display for BlockTensor
function Base.show(io::IO, ::MIME"text/plain", bt::BlockTensor)
    println(io, "BlockTensor with group: ", bt.group)
    println()
    println(io, "Dimensions: ")
    for (i, obj_dict) in enumerate(bt.objects)
        println(io, "  Leg $i: $obj_dict")
    end
    println()
    println(io, "Blocks: ")
    println()
    for (key, value) in bt.data
        println(io, "  Sector: ", key, "\n  Block:")
        Base.show(io, MIME("text/plain"), value)
        println()
        println()
    end
end

# Function to set a specific block of the tensor
function set_block!(bt::BlockTensor{Float64, G}, sector::Vector{GroupElement{G}}, value::Array{Float64}) where {G <: Group}
    # Ensure the key matches the group and dimensions are consistent
    if length(sector) != length(bt.objects)
        throw(ArgumentError("The sector length does not match the number of tensor legs."))
    end

    if multiply(sector) != identity(bt.group)
        throw(ArgumentError("The sector multiplies to non identity."))
    end

    for (i, g) in enumerate(sector)
        if bt.objects[i].obj[g] == 0
            throw(ArgumentError("Group element $g is not valid for leg $i."))
        end
        expected_dim = bt.objects[i].obj[g]
        if size(value, i) != expected_dim
            throw(ArgumentError("The size of the array along dimension $i does not match the expected dimension $expected_dim for group element $g."))
        end
    end
    bt.data[sector] = value
    return
end

# D4 = DihedralGroup(4)
# e = GroupElement((0, 0), D4)
# s = GroupElement((1, 0), D4)
# r = GroupElement((0,1),D4)
# r2 = GroupElement((0, 2), D4)
# sr2 = s * r2
# A = Object(Dict(e => 2, s*r => 3, s => 1, r2 => 2, sr2 => 1, r => 5), D4)

# @show T = BlockTensor(D4, [A, A, A, A])

# set_block!(T, [e, s*r, s, r], rand(2, 3, 1, 5))
# set_block!(T, [r2, r2, s, s], rand(2,2,1,1))
# T

# x = 123
# # Function to get a block of the tensor
# function get_block(bt::BlockTensor{T, G}, key::Vector{G}) where {T, G <: GroupElement}
#     return get(bt.data, key, nothing)
# end





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


include("groups.jl")


struct Object{G <: Group, T}
    obj::Dict{GroupElement{G, T}, Int}
    group::G
    function Object(X::Dict{GroupElement{G, T}, Int}, group::G) where {G<:Group, T}
        el = elements(group)
        for g in el
            if !haskey(X, g)
                X[g] = 0
            end
        end
        return new{G, T}(X, group)
    end
end

# Custom display for Object
function Base.show(io::IO, x::Object{G}) where {G<:Group}
    group = x.group
    obj = x.obj
    el = elements(group)
    j = 1
    for (j, i) in enumerate(el)
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

# Custom display for a sector
function Base.show(io::IO, v::Tuple{Vararg{GroupElement{G}}}) where {G<:Group}
    for (i, g) in enumerate(v)
        if i != 1
            print(io, " ⊗ ")
        end
        print(io, g)
    end
end

# Define BlockTensor struct (generic value type T)
struct BlockTensor{T, G <: Group}
    group::G
    objects::Vector{Object{G}}
    data::Dict{Tuple{Vararg{GroupElement{G}}}, Array{T}}
end

# Constructor
function BlockTensor(group::G, objects::Vector{Object{G}}) where {G<:Group}
    data = Dict{Tuple{Vararg{GroupElement{G}}}, Array{Any}}()
    return BlockTensor{Any, G}(group, objects, data)
end

# Custom display for BlockTensor
function Base.show(io::IO, ::MIME"text/plain", bt::BlockTensor)
    println(io, "BlockTensor with group: ", bt.group)
    println(io, "Number of legs: ", length(bt.objects))
    println(io, "Number of blocks: ", length(bt.data))
    println("\nBlocks: ")
    max_display = 5
    for (i, (key, value)) in enumerate(bt.data)
        if i > max_display
            println(io, "... (and $(length(bt.data) - max_display) more blocks)")
            break
        end
        println(io, "  Sector: ", key, " → Block size: ", size(value))
    end
end

# Set a specific block
function set_block!(bt::BlockTensor{T, G}, sector::Tuple{Vararg{GroupElement{G}}}, value::Array{T}) where {T, G <: Group}
    if length(sector) != length(bt.objects)
        throw(ArgumentError("Sector length $(length(sector)) does not match number of tensor legs $(length(bt.objects))."))
    end
    if multiply(sector) != identity(bt.group)
        throw(ArgumentError("Sector $sector does not multiply to the identity."))
    end
    for (i, g) in enumerate(sector)
        if bt.objects[i].obj[g] == 0
            throw(ArgumentError("Group element $g is not valid for leg $i. Valid elements: $(keys(bt.objects[i].obj))"))
        end
        expected_dim = bt.objects[i].obj[g]
        if size(value, i) != expected_dim
            throw(ArgumentError("Dimension mismatch for leg $i with group element $g. Expected: $expected_dim, got: $(size(value, i))"))
        end
    end
    bt.data[sector] = value
end

# # Get a block
# function get_block(bt::BlockTensor{T, G}, sector::Tuple{Vararg{GroupElement{G}}}) where {T, G <: Group}
#     return get(bt.data, sector, nothing)
# end

# # Add two BlockTensors (block-wise addition)
# function +(bt1::BlockTensor{T, G}, bt2::BlockTensor{T, G}) where {T, G <: Group}
#     if bt1.group != bt2.group || bt1.objects != bt2.objects
#         throw(ArgumentError("BlockTensors must have the same group and objects to add them."))
#     end
#     result = BlockTensor(bt1.group, bt1.objects)
#     for (key, value) in bt1.data
#         result.data[key] = copy(value)
#     end
#     for (key, value) in bt2.data
#         if haskey(result.data, key)
#             result.data[key] .+= value
#         else
#             result.data[key] = copy(value)
#         end
#     end
#     return result
# end





include("test.jl")
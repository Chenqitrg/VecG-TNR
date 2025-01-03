include("groups.jl")


struct Obj
    sumd::Dict{GroupElement{<:Group}, Int}  # 直接使用 GroupElement 的子类型作为键的类型
    function Obj(pairs::Pair{GroupElement{<:Group}, Int}...)
        g0 = first(pairs).first
        group = g0.group
        sumd = Dict{GroupElement{<:Group}, Int}()
        for pair in pairs
            sumd[pair.first] = pair.second
        end
        el = elements(group)
        for g in el
            if !haskey(sumd, g)
                sumd[g] = 0
            end
        end
        return new(sumd)
    end
end

# # Custom display for Object
# function Base.show(io::IO, x::Object{G}) where {G<:Group}
#     group = x.group
#     obj = x.obj
#     el = elements(group)
#     j = 1
#     for (j, i) in enumerate(el)
#         if obj[i] != 0
#             if j != 1
#                 print(io, " ⊕ ")
#             end
#             if obj[i] == 1
#                 print(io, i)
#             else
#                 print(io, obj[i], i)
#             end
#         end
#     end
# end

# # Custom display for a sector
# function Base.show(io::IO, v::Tuple{Vararg{GroupElement{G}}}) where {G<:Group}
#     for (i, g) in enumerate(v)
#         if i != 1
#             print(io, " ⊗ ")
#         end
#         print(io, g)
#     end
# end

# # Define BlockTensor struct (generic value type T)
# struct BlockTensor{T, G <: Group}
#     group::G
#     objects::Tuple{Vararg{Object{G}}}
#     data::Dict{Tuple{Vararg{GroupElement{G}}}, Array{T}}
# end

# # Constructor
# function BlockTensor(element_type::Type, group::G, objects::Tuple{Vararg{Object{G}}}) where {G<:Group}
#     data = Dict{Tuple{Vararg{GroupElement{G}}}, Array{Any}}()
#     return BlockTensor{element_type, G}(group, objects, data)
# end

# # Custom display for BlockTensor
# function Base.show(io::IO, ::MIME"text/plain", bt::BlockTensor)
#     println(io, "BlockTensor with group: ", bt.group)
#     println(io, "Number of legs: ", length(bt.objects))
#     for (i, obj) in enumerate(bt.objects)
#         println(io, "Leg $i is of object $obj")
#     end
#     println("Blocks: ")
#     max_display = 5
#     for (i, (key, value)) in enumerate(bt.data)
#         if i > max_display
#             println(io, "... (and $(length(bt.data) - max_display) more blocks)")
#             break
#         end
#         println(io, "  Sector: ", key, " → Block size: ", size(value))
#     end
# end

# # Set a specific block
# function set_block!(bt::BlockTensor{T, G}, sector::Tuple{Vararg{GroupElement{G}}}, value::Array{T}) where {T, G <: Group}
#     if length(sector) != length(bt.objects)
#         throw(ArgumentError("Sector length $(length(sector)) does not match number of tensor legs $(length(bt.objects))."))
#     end
#     if multiply(sector) != identity(bt.group)
#         throw(ArgumentError("Sector $sector does not multiply to the identity."))
#     end
#     for (i, g) in enumerate(sector)
#         if bt.objects[i].obj[g] == 0
#             throw(ArgumentError("Group element $g is not valid for leg $i."))
#         end
#         expected_dim = bt.objects[i].obj[g]
#         if size(value, i) != expected_dim
#             throw(ArgumentError("Dimension mismatch for leg $i with group element $g. Expected: $expected_dim, got: $(size(value, i))"))
#         end
#     end
#     bt.data[sector] = value
# end

# function sector_size(bt::BlockTensor{T, G}, sector::Tuple{Vararg{GroupElement{G}}}) where {T, G <: Group}
#     group = bt.group
#     if multiply(sector) != identity(group)
#         throw(ArgumentError("The sector $sector is not consistent"))
#     else
#         size = ()
#         for (i, g) in enumerate(sector)
#             object = bt.objects[i]
#             multiplicity = object.obj[g]
#             size = (size..., multiplicity)
#         end
#         return size
#     end
# end

# function random_block_tensor(element_type::Type, group::G, objects::Tuple{Vararg{Object{G}}}) where {G<:Group}
#     bt = BlockTensor(element_type, group, objects)
#     e = identity(group)
#     iter = group_tree(group, e, length(objects))
#     for sector in iter
#         size = sector_size(bt, sector)
#         bt.data[sector] = rand(element_type, size...)
#     end
#     return bt
# end





# # include("test.jl")
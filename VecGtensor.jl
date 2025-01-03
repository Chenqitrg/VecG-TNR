include("groups.jl")

abstract type UFC end
abstract type PointedCat <: UFC end

struct VecG{G<:Group}
    group::G
end


struct Obj{G<:Group}
    sumd::Dict{GroupElement{G}, Int}  # 直接使用 GroupElement 的子类型作为键的类型
end

function Obj(pairs::Pair{GroupElement{G}, Int}...) where G<:Group
    g0 = first(pairs).first
    group = g0.group
    sumd = Dict{GroupElement{G}, Int}()
    for pair in pairs
        sumd[pair.first] = pair.second
    end
    el = elements(group)
    for g in el
        if !haskey(sumd, g)
            sumd[g] = 0
        end
    end
    return Obj(sumd)
end

function zero_obj(group::G) where G<:Group
    sumd = Dict{GroupElement{G}, Int}()
    el = elements(group)
    for g in el
        sumd[g] = 0
    end
    return Obj(sumd)
end

function get_group(ob::Obj)
    return first(keys(ob.sumd)).group
end


# Custom display for Object
function Base.show(io::IO, x::Obj)
    group = get_group(x)
    sumd = x.sumd
    el = elements(group)
    op = false
    for i in el
        if sumd[i] != 0
            if op == true
                print(io, " ⊕ ")
            end
            if sumd[i] == 1
                print(io, i)
            else
                print(io, sumd[i], i)
            end
            op = true
        end
    end
end

function Base.getindex(obj::Obj{G}, key::GroupElement{G}) where G<:Group
    return get(obj.sumd, key, 0)  # 如果键不存在，默认返回 0
end


struct Sector{G<:Group}
    sect::Tuple{Vararg{GroupElement{G}}}
end

function Sector(groupelements::GroupElement...)
    group = first(groupelements).group
    if multiply(groupelements) != identity(group)
        throw(ArgumentError("The sector $groupelements is not consistent"))
    else
        return Sector(groupelements)
    end
end


function Base.getindex(S::Sector, key::Int)
    return get(S.sect, key, 0)  # 如果键不存在，默认返回 0
end


# Custom display for a sector
function Base.show(io::IO, S::Sector)
    for (i, g) in enumerate(S.sect)
        if i != 1
            print(io, " ⊗ ")
        end
        print(io, g)
    end
end

D4 = DihedralGroup(4)
s = GroupElement((1,0),D4)
r = GroupElement((0,1),D4)
e = identity(D4)
A = Obj(e=>1, s=>2, r=>3)
get_group(A)
zero_obj(D4)

P = Sector(e,s,r)

Z4 = CyclicGroup(4)
a = GroupElement(1,Z4)
@show S = Sector(a, a*a)

# Define BlockTensor struct (generic value type T)
struct Mor{G <: Group, T}
    objects::Tuple{Vararg{Obj{G}}}
    data::Dict{Sector{G}, Array{T}}
end

# Constructor
function Mor(element_type::Type, objects::Tuple{Vararg{Obj{G}}}) where {G<:Group}
    data = Dict{Tuple{Vararg{GroupElement{G}}}, Array{Any}}()
    return Mor{G,element_type}(objects, data)
end

function get_group(T::Mor)
    return get_group(T.objects[1])
end

function Base.getindex(T::Mor, g::GroupElement...)
    S = Sector(g)
    return T.data[S]
end

function Base.getindex(T::Mor, i::Int)
    return T.objects[i]
end

T = Mor(Float64, (A,A,A))

# Custom display for Morphism
function Base.show(io::IO, ::MIME"text/plain", T::Mor)
    println(io, "Group: ", get_group(T))
    for (i, obj) in enumerate(T.objects)
        println(io, "Leg $i is of object $obj")
    end
end


# Set a specific block
function Base.setindex!(mor::Mor{G, T}, value::Array{T}, g::GroupElement{G}...) where {T, G <: Group}
    terminal_size = sector_size(mor, g)
    data_size = size(value)
    if terminal_size != data_size
        throw(ArgumentError("Data size $data_size does not match sector size $terminal_size"))
    end
    mor.data[Sector(g)] = value
end

function sector_size(mor::Mor{G, T}, tup::Tuple{Vararg{GroupElement{G}}}) where {T, G<:Group}
    group = get_group(mor)
    if multiply(tup) != identity(group)
        throw(ArgumentError("The sector $tup is not consistent"))
    else
        size = ()
        for (i, g) in enumerate(tup)
            object = mor[i]
            multiplicity = object[g]
            size = (size..., multiplicity)
        end
        return size
    end
end

function sector_size(mor::Mor{G, T}, sector::Sector{G}) where {T, G<:Group}
    return sector_size(mor, sector.sect)
end

function random_mor(element_type::Type, objects::Tuple{Vararg{Obj}})
    mor = Mor(element_type, objects)
    group = get_group(mor)
    e = identity(group)
    iter = group_tree(e, length(objects))
    for grouptuple in iter
        size = sector_size(mor, grouptuple)
        mor[grouptuple...] = rand(element_type, size...)
    end
    return mor
end

T = random_mor(Float64, (A, A, A))



# # include("test.jl")
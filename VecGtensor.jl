# include("groups.jl")
# include("display.jl")
abstract type UFC end
abstract type PointedCat <: UFC end

struct VecG{G<:Group}
    group::G
end

struct Obj{G<:Group}
    sumd::Dict{GroupElement{G}, Int}  # 直接使用 GroupElement 的子类型作为键的类型
end

struct Sector{G<:Group}
    sect::Tuple{Vararg{GroupElement{G}}}
end

struct Mor{G <: Group, T}
    objects::Tuple{Vararg{Obj{G}}}
    data::Dict{Sector{G}, Array{T}}
end

# Construction function of an object in VecG
function Obj(pairs::Pair{GroupElement{G}, Int}...) where G <: Group
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

# Construction function for a sector
function Sector(groupelements::GroupElement...)
    group = first(groupelements).group
    if multiply(groupelements) != identity(group)
        throw(ArgumentError("The sector $groupelements is not consistent"))
    else
        return Sector(groupelements)
    end
end

# Construction function for a morphism
function Mor(element_type::Type, objects::Tuple{Vararg{Obj{G}}}) where {G<:Group}
    data = Dict{Tuple{Vararg{GroupElement{G}}}, Array{Any}}()
    return Mor{G,element_type}(objects, data)
end

# Construct a zero object
function zero_obj(group::G) where G<:Group
    sumd = Dict{GroupElement{G}, Int}()
    el = elements(group)
    for g in el
        sumd[g] = 0
    end
    return Obj(sumd)
end

function dual_obj(obj::Obj{G}) where G<:Group
    dualobj = Dict{GroupElement{G}, Int}()
    dict = obj.sumd
    for g in keys(dict)
        dualobj[g] = dict[inverse(g)]
    end
    return Obj(dualobj)
end

function Base.:(==)(X::Obj{G}, Y::Obj{G}) where G<:Group
    return X.sumd == Y.sumd
end

# Construct a random morphism for a given objects
function random_mor(element_type::Type, objects::Tuple{Vararg{Obj}})
    mor = Mor(element_type, objects)
    group = get_group(mor)
    e = identity(group)
    iter = group_tree(e, length(objects))
    for grouptuple in iter
        size = get_sector_size(mor, grouptuple)
        mor[grouptuple...] = rand(element_type, size...)
        # @show grouptuple
        # @show size
    end
    return mor
end

# Get the group that the object belongs to
function get_group(ob::Obj)
    return first(keys(ob.sumd)).group
end

# Get the direct summand multiplicity
function Base.getindex(obj::Obj{G}, key::GroupElement{G}) where G<:Group
    return get(obj.sumd, key, 0)  # 如果键不存在，默认返回 0
end

function Base.setindex!(obj::Obj{G}, multiplicity::Int, key::GroupElement{G}) where G<:Group
    obj.sumd[key] = multiplicity
end

# Get the group of a morphism
function get_group(T::Mor)
    return get_group(T.objects[1])
end

# Get the key-th group element
function Base.getindex(S::Sector, key::Int)
    return get(S.sect, key, 0)  # 如果键不存在，默认返回 0
end

# Get the object of the i-th leg
function Base.getindex(T::Mor, i::Int)
    return T.objects[i]
end

function Base.lastindex(T::Mor)
    return length(T.objects)
end


# Get the object of the leg within the range
function Base.getindex(mor::Mor, range::UnitRange{Int})
    return mor.objects[range]
end

# Get the tensor of the g... sector
function Base.getindex(mor::Mor, g::GroupElement...)
    S = Sector(g)
    return mor.data[S]
end

# Get the tensor of the S sector
function Base.getindex(mor::Mor{G, T}, S::Sector{G}) where {T, G<:Group}
    return mor.data[S]
end

# Get the size of a morphism at a given sector
function get_sector_size(mor::Mor{G, T}, tup::Tuple{Vararg{GroupElement{G}}}) where {T, G<:Group}
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

# Get the size of the tensor for a given sector
function get_sector_size(mor::Mor{G, T}, sector::Sector{G}) where {T, G<:Group}
    return get_sector_size(mor, sector.sect)
end

function Base.setindex!(mor::Mor{G, T}, obj::Obj{G}, i::Int) where {T, G<:Group}
    group = get_group(obj)
    for g in elements(group)
        mor[i][g] = obj[g]
    end
end


# Set a specific block
function Base.setindex!(mor::Mor{G, T}, value::Array{T}, g::GroupElement{G}...) where {T, G <: Group}
    terminal_size = get_sector_size(mor, g)
    data_size = size(value)
    if terminal_size != data_size
        throw(ArgumentError("Data size $data_size does not match sector size $terminal_size"))
    end
    mor.data[Sector(g)] = value
end

# Set a specific sector
function Base.setindex!(mor::Mor{G, T}, value::Array{T}, S::Sector{G}) where {T, G <: Group}
    terminal_size = get_sector_size(mor, S)
    data_size = size(value)
    if terminal_size != data_size
        throw(ArgumentError("Data size $data_size does not match sector size $terminal_size"))
    end
    mor.data[S] = value
end

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
    if multiply(groupelements) != identity_element(group)
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
    e = identity_element(group)
    iter = group_tree(e, length(objects))
    for grouptuple in iter
        size = get_sector_size(mor, grouptuple)
        mor[grouptuple...] = rand(element_type, size...)
    end
    return mor
end

# Construct a random morphism for a given objects
function zero_mor(element_type::Type, objects::Tuple{Vararg{Obj}})
    mor = Mor(element_type, objects)
    group = get_group(mor)
    e = identity_element(group)
    iter = group_tree(e, length(objects))
    for grouptuple in iter
        size = get_sector_size(mor, grouptuple)
        mor[grouptuple...] = zeros(element_type, size...)
    end
    return mor
end

function identity_mor(element_type::Type, object::Obj)
    dual = dual_obj(object)
    group = get_group(object)
    delta = Mor(element_type, (object, dual))
    for g in elements(group)
        sect = get_sector_size(delta, (g, inverse(g)))
        delta[g, inverse(g)] = Matrix{element_type}(I, sect...)
    end

    return delta
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
    if multiply(tup) != identity_element(group)
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


function is_accend(tup::Tuple{Vararg{Int}}, modn::Int)
    test = true
    for i = 1 : length(tup)-1
        if !(tup[i] in 1:modn) || !(tup[i+1] in 1:modn) || tup[i+1] != (mod(tup[i], modn) + 1)
            test = false
        end
    end
    return test
end

function is_decend(tup::Tuple{Vararg{Int}}, modn::Int)
    test = true
    for i = 1 : length(tup)-1
        if !(tup[i] in 1:modn) || !(tup[i+1] in 1:modn) || tup[i+1] != (mod(tup[i], modn) - 1)
            test = false
        end
    end
    return test
end

function is_cyclic(tup::Tuple{Vararg{Int}}, modn::Int)
    test = is_accend(tup, modn)
    if tup[1]!=(mod(tup[end], modn) + 1)
        test = false
    end
    return test
end

function to_perm(tup::Tuple{Vararg{Int}}, modn::Int)
    if !is_accend(tup, modn)
        throw(Argumenterror("The $tup is not in an accending order"))
    end
    newtup = tup[1]:tup[1]+modn-1
    mod_newtup = Tuple(map(x->mod(x-1,modn)+1, newtup))
    return mod_newtup
end

function VecG_permutesectors(sect::Sector, perm::Tuple{Vararg{Int}})
    modn = length(sect.sect)
    if is_cyclic(perm, modn) == false
        throw(ArgumentError("The permutation $perm is not cyclic"))
    end
    perm_sect_tup = ()
    for i in perm
        perm_sect_tup = (perm_sect_tup..., sect[i])
    end
    return Sector(perm_sect_tup...)
end

function VecG_permutedims(mor::Mor{G, T}, perm::Tuple{Vararg{Int}}) where {T, G<:Group}
    modn = length(mor.objects)
    if is_cyclic(perm, modn) == false
        throw(ArgumentError("The permutation $perm is not cyclic"))
    end
    perm_obj = ()
    for perm_i in perm
        perm_obj = (perm_obj..., mor[perm_i])
    end
    perm_mor = Mor(T, perm_obj)
    for sect in keys(mor.data)
        perm_sect = VecG_permutesectors(sect, perm)
        perm_tensor = permutedims(mor[sect], perm)
        perm_mor[perm_sect] = perm_tensor
    end
    return perm_mor
end

function Base.broadcasted(::typeof(/), mor::Mor, x::Number)
    for key in keys(mor.data)
        mor[key] = mor[key] ./ x
    end
    return mor
end

function Base.broadcasted(::typeof(/), x::Number, mor::Mor{G, T}) where {T, G<:Group}
    newmor = Mor(T, mor.objects)
    for key in keys(mor.data)
        newmor[key] = diagm(x ./ diag(mor[key]))
    end
    return newmor
end

function Base.broadcasted(::typeof(sqrt), mor::Mor{G, T}) where {T, G<:Group}
    newmor = Mor(T, mor.objects)
    for key in keys(mor.data)
        newmor[key] = sqrt.(mor[key])
    end
    return newmor
end

#   | | | | |
#   1 2 3 4 5
#   | | | | |
#   ^ ^ ^ ^ ^
#   | | | | |
#       T
# dagger:
#   | | | | |
#   5 4 3 2 1
#   | | | | |
#   v v v v v
#   | | | | |
#       T^†
function VecG_dag(mor::Mor{G, T}) where {T, G<:Group}
    group = get_group(mor)
    leg_number = length(mor.objects)
    obj = ()
    for i in 1:leg_number
        obj = (dual_obj(mor[i]), obj...)
    end
    mor_dag = Mor(T, obj)
    for key in group_tree(identity_element(group), leg_number)
        # @show key
        dag_key = reverse(inverse.(key))
        mor_dag[key...] = permutedims(conj.(mor[dag_key...]), reverse(1:leg_number))
    end
    return mor_dag
end

function max_abs(mor::Mor)
    max = 0.
    for key in keys(mor.data)
        if !(0 in size(mor[key]))
            localmax = maximum(abs.(mor[key]))
            if localmax > max
                max = localmax
            end
        end
    end
    return max
end

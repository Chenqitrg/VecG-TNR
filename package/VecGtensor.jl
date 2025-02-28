abstract type UFC end
abstract type PointedCat <: UFC end

struct VecG{G<:Group}
    group::G
end


"""
Structure of objects in VecG. To be used in Mor.

To adapt to infinite groups, we use a dictionary to store the multiplicities of group elements. Moreover, some group elements may not appear in the object.

# Input: 
- pairs of the form g=>n

# Output:
- an object of the form n1g1⊕n2g2⊕n3g3...

# Example

```
julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       Obj(e=>1, a=>2, aa=>3)

e⊕2a⊕3a²
```

"""
struct Obj{G<:Group}
    sumd::Dict{GroupElement{G}, Int}  # 直接使用 GroupElement 的子类型作为键的类型
end
function Obj(pairs::Pair{GroupElement{G}, Int}...) where G <: Group
    sumd = Dict{GroupElement{G}, Int}()
    for pair in pairs
        if pair.second == 0
            continue
        else
            sumd[pair.first] = pair.second
        end
    end
    return Obj(sumd)
end

"""
Generating a zero object.

# Input:
- a group

# Output:
- an object with all zero multiplicities

# Example

```
julia> G = CyclicGroup(3)
       zero_obj(G)

"""
function zero_obj(group::G) where G<:Group
    sumd = Dict{GroupElement{G}, Int}()
    el = identity_element(group)
    sumd[el] = 0
    return Obj(sumd)
end

"""
Generating a dual object.

# Input:
- an object

# Output:
- the dual object, whose multiplicities are the same as the input object, but the group elements are inversed.

# Example

```
julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       A = Obj(e=>1, a=>2, aa=>3)
       dual_obj(A)

e⊕3a⊕2a²

```
"""
function dual_obj(obj::Obj{G}) where G<:Group
    dualobj = Dict{GroupElement{G}, Int}()
    dict = obj.sumd
    for g in keys(dict)
        dualobj[inverse(g)] = dict[g]
    end
    return Obj(dualobj)
end

"""
Comparing two objects.
"""
function Base.:(==)(X::Obj{G}, Y::Obj{G}) where G<:Group
    return X.sumd == Y.sumd
end

"""
Structure of sectors in VecG

# Input:
- a tuple of group elements

# Output:
- a sector of the form g1⊗g2⊗g3..., such that g1g2g3... = e
If the sector is not consistent, an error will be thrown.

# Example

```
julia> G = CyclicGroup(3)
       e = GroupElement((0,0), G)
       r = GroupElement((0,1), G)
       s = GroupElement((1,0), G)
       Sector(s, r, s*r, e)

s⊗r⊗sr⊗e

```
"""
struct Sector{G<:Group}
    sect::Tuple{Vararg{GroupElement{G}}}
end
function Sector(groupelements::GroupElement...)
    group = first(groupelements).group
    if multiply(groupelements) != identity_element(group)
        throw(ArgumentError("The sector $groupelements is not consistent"))
    else
        return Sector(groupelements)
    end
end

"""
The key data type in VecG_TNR.

# Components:
- objects: a tuple of objects
- data: a dictionary, whose keys are sectors and values are tensors

# Construction function:
- Mor(element_type, objects)
## Input:
- element type: data type, for example Float64, ComplexF64, etc...
- objects: a tuple of objects

## Output: 
- an empty tensor, stored by a dictionary, whose key are sectors

"""
struct Mor{G <: Group, T}
    objects::Tuple{Vararg{Obj{G}}}
    data::Dict{Sector{G}, Array{T}}
end
function Mor(element_type::Type, objects::Tuple{Vararg{Obj{G}}}) where {G<:Group}
    data = Dict{Sector{G}, Array{element_type}}() # Previously the type of the data here does not match the type of data in Mor, but I do not know why there is no error
    return Mor{G,element_type}(objects, data)
end

"""
Generating a random morphism for a given objects.

# Input:
- element type: data type, for example Float64, ComplexF64, etc...
- objects: a tuple of objects

# Output:
- a random morphism, whose data are random numbers

# Example

```
julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       A = Obj(e=>1, a=>2, aa=>3)
       B = Obj(e=>2, a=>3, aa=>2)
       C = Obj(e=>1, a=>2, aa=>3)
       D = Obj(e=>2, a=>3, aa=>2)
       T = random_mor(Float64, (A, B, C, D))
         The data is too long to show. 
```
"""
function random_mor(element_type::Type, objects::Tuple{Vararg{Obj}})
    mor = Mor(element_type, objects)
    iter = Iterators.product(map(x->keys(x.sumd), objects)...)
    group = get_group(mor)
    e = identity_element(group)
    iter = group_tree(e, length(objects))
    for grouptuple in iter
        size = get_sector_size(mor, grouptuple)
        mor[grouptuple...] = rand(element_type, size...)
    end
    return mor
end

"""
Generating a zero morphism for a given objects.

# Input:
- element type: data type, for example Float64, ComplexF64, etc...
- objects: a tuple of objects

# Output:
- a zero morphism, whose data are all zeros

# Example

```
julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       A = Obj(e=>1, a=>2, aa=>3)
       B = Obj(e=>2, a=>3, aa=>2)
       C = Obj(e=>1, a=>2, aa=>3)
       D = Obj(e=>2, a=>3, aa=>2)
       T = zero_mor(Float64, (A, B, C, D))
         The data is too long to show.
```
"""
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

"""
Generating an identity morphism for a given object.

# Input:
- element type: data type, for example Float64, ComplexF64, etc...
- object: an object

# Output:
- an identity morphism, whose data are identity matrices

# Example

```
julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       A = Obj(e=>1, a=>2, aa=>3)
       identity_mor(Float64, A)
```
"""
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
"""
Get the group that the object belongs to

# Input:
- an object

# Output:
- the group that the object belongs to

# Example

```
julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       A = Obj(e=>1, a=>2, aa=>3)
       get_group(A)
```
"""
function get_group(ob::Obj)
    return first(keys(ob.sumd)).group
end

"""
Get the multiplicity of a group element in an object

# Input:
- an object
- a group element

# Output:
- the multiplicity of the group element in the object

# Example

```
julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       A = Obj(e=>1, a=>2, aa=>3)
       A[a]
       2
```
"""
function Base.getindex(obj::Obj{G}, key::GroupElement{G}) where G<:Group
    return get(obj.sumd, key, 0)  # 如果键不存在，默认返回 0
end


"""
Set the multiplicity of a group element in an object

# Input:
- an object
- the multiplicity
- a group element

# Output:
- the object with the multiplicity of the group element set

# Example

```

julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       A = Obj(e=>1, a=>2, aa=>3)
       A[a] = 5
       e⊕5a⊕3a²
```
"""
function Base.setindex!(obj::Obj{G}, multiplicity::Int, key::GroupElement{G}) where G<:Group
    obj.sumd[key] = multiplicity
end


"""
Get the group of a morphism

# Input:
- a morphism

# Output:
- the group of the morphism

# Example

```

julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       A = Obj(e=>1, a=>2, aa=>3)
       B = Obj(e=>2, a=>3, aa=>2)
       C = Obj(e=>1, a=>2, aa=>3)
       D = Obj(e=>2, a=>3, aa=>2)
       T = Mor(Float64, (A, B, C, D))
       get_group(T)
       ℤ₃
```
"""
function get_group(T::Mor)
    return get_group(T.objects[1])
end


"""
Get the key-th group element in a sector

# Input:
- a sector
- the index

# Output:
- the key-th group element in the sector

# Example

```

julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       S = Sector(e, a, aa)
       S[2]
       a
```
"""
function Base.getindex(S::Sector, key::Int)
    return get(S.sect, key, 0)  # 如果键不存在，默认返回 0
end

# Get the object of the i-th leg
"""
Get the object of the i-th leg

# Input:
- a morphism
- the index

# Output:
- the object of the i-th leg

# Example

```
julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       A = Obj(e=>1, a=>2, aa=>3)
       B = Obj(e=>2, a=>3, aa=>2)
       C = Obj(e=>1, a=>2, aa=>3)
       D = Obj(e=>2, a=>3, aa=>2)
       T = Mor(Float64, (A, B, C, D))
       T[1]
       e ⊕ 2a ⊕ 3a²
```
"""
function Base.getindex(T::Mor, i::Int)
    return T.objects[i]
end

"""
Get the number of legs in a morphism

# Input:
- a morphism

# Output:
- the number of legs in the morphism

# Example

```
julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       A = Obj(e=>1, a=>2, aa=>3)
       B = Obj(e=>2, a=>3, aa=>2)
       C = Obj(e=>1, a=>2, aa=>3)
       D = Obj(e=>2, a=>3, aa=>2)
       T = Mor(Float64, (A, B, C, D))
        lastindex(T)
       4
```
"""
function Base.lastindex(T::Mor)
    return length(T.objects)
end


"""
Get a vector of objects in a morphism

# Input:
- a morphism
- a range

# Output:
- a vector of objects in the morphism

# Example

```
julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       A = Obj(e=>1, a=>2, aa=>3)
       B = Obj(e=>2, a=>3, aa=>2)
       C = Obj(e=>1, a=>2, aa=>3)
       D = Obj(e=>2, a=>3, aa=>2)
       T = Mor(Float64, (A, B, C, D))
       T[1:2]
       (e ⊕ 2a ⊕ 3a², e ⊕ 3a ⊕ 2a²)
```
"""
function Base.getindex(mor::Mor, range::UnitRange{Int})
    return mor.objects[range]
end

"""
Get the tensor of the sector

# Input:
- a morphism
- a sector

# Output:
- the tensor of the sector

Both the form T[S] and T[S.sect], which can be, for example T[g,h,k,l], are supported.

# Example

```
julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       A = Obj(e=>1, a=>2, aa=>3)
       B = Obj(e=>2, a=>3, aa=>2)
       C = Obj(e=>1, a=>2, aa=>3)
       D = Obj(e=>2, a=>3, aa=>2)
       T = Mor(Float64, (A, B, C, D))
       S = Sector(e, a, aa)
       T[S]
       T[e,a,aa]
```

"""
function Base.getindex(mor::Mor, g::GroupElement...)
    S = Sector(g)
    return mor.data[S]
end
function Base.getindex(mor::Mor{G, T}, S::Sector{G}) where {T, G<:Group}
    return mor.data[S]
end

"""
Get the size of the tensor for a given sector

# Input:
- a morphism
- a tuple of group elements or a sector

# Output:
- the size of the tensor for the sector constructed by the tuple of group elements

# Example

```
julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       A = Obj(e=>1, a=>2, aa=>3)
       B = Obj(e=>2, a=>3, aa=>2)
       C = Obj(e=>1, a=>2, aa=>3)
       D = Obj(e=>2, a=>3, aa=>2)
       T = random_mor(Float64, (A, B, C, D))
       get_sector_size(T, (e, a, aa, e))
       (1, 3, 3, 2)

julia> S = Sector(e, a, aa, e)
    get_sector_size(T, S)
    (1, 3, 3, 2)
```
"""
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
function get_sector_size(mor::Mor{G, T}, sector::Sector{G}) where {T, G<:Group}
    return get_sector_size(mor, sector.sect)
end

# So far it seems to be a good idea to delete the following function
# function Base.setindex!(mor::Mor{G, T}, obj::Obj{G}, i::Int) where {T, G<:Group}
#     group = get_group(obj)
#     for g in elements(group)
#         mor[i][g] = obj[g]
#     end
# end


"""
Set the morphism to a given array for a given sector

# Input:
- a morphism
- an array
- a tuple of group elements or a sector

# Output:
- the morphism with the tensor set to the array for the sector constructed by the tuple of group elements

# Example

```
julia> G = CyclicGroup(3)
       e = GroupElement(0, G)
       a = GroupElement(1, G)
       aa = GroupElement(2, G)
       A = Obj(e=>1, a=>2, aa=>3)
       B = Obj(e=>2, a=>3, aa=>2)
       C = Obj(e=>1, a=>2, aa=>3)
       D = Obj(e=>2, a=>3, aa=>2)
       T = random_mor(Float64, (A, B, C, D))
       T[e,a,aa,e] = rand(1,3,3,2)
julia> S = Sector(e, a, aa, e)
       T[S] = rand(1,3,3,2)
```
"""
function Base.setindex!(mor::Mor{G, T}, value::Array{T}, g::GroupElement{G}...) where {T, G <: Group}
    terminal_size = get_sector_size(mor, g)
    data_size = size(value)
    if terminal_size != data_size
        throw(ArgumentError("Data size $data_size does not match sector size $terminal_size"))
    end
    mor.data[Sector(g)] = value
end
function Base.setindex!(mor::Mor{G, T}, value::Array{T}, S::Sector{G}) where {T, G <: Group}
    terminal_size = get_sector_size(mor, S)
    data_size = size(value)
    if terminal_size != data_size
        throw(ArgumentError("Data size $data_size does not match sector size $terminal_size"))
    end
    mor.data[S] = value
end

# Here is bug
"""
Judge whether a tuple is in an accending order.

# Input:
- a tuple of integers
- a modulus

# Output:
- a boolean value, indicating whether the tuple is in an accending order
It reports false if the tuple is not in an accending order or the elements are not in the range of 1 to modn.

# Example

```
julia> is_accend((1,2,3,4), 4)
       true

julia> is_accend((4,1), 4)
    true
```
"""
function is_accend(tup::Tuple{Vararg{Int}}, modn::Int)
    test = true
    for i = 1 : length(tup)-1
        if !(tup[i] in 1:modn) || !(tup[i+1] in 1:modn) || tup[i+1] != (mod(tup[i], modn) + 1)
            test = false
        end
    end
    return test
end

"""
Judge whether a tuple is in a descending order.

# Input:
- a tuple of integers
- a modulus

# Output:
- a boolean value, indicating whether the tuple is in a descending order
It reports false if the tuple is not in a descending order or the elements are not in the range of 1 to modn.

# Example

```
julia> is_descend((4,3,2,1), 4)
       true

julia> is_descend((4,1), 4) 
         false

```
"""
function is_descend(tup::Tuple{Vararg{Int}}, modn::Int)
    test = true
    for i = 1 : length(tup)-1
        if (!(tup[i] in 1:modn)) || (!(tup[i+1] in 1:modn)) || (mod(tup[i+1], modn) != (tup[i] - 1))
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

function Base.:+(mora::Mor{G, T}, morb::Mor{G, T}) where {T, G<:Group}
    if length(mora.objects)!=length(morb.objects)
        throw(ArgumentError("The first tensor has $(length(mora.objects)) number of legs, while the second tensor has $(length(morb.objects)) number of legs"))
    else
        for i in 1:length(mora.objects)
            if mora[i]!=morb[i]
                throw(ArgumentError("The $i-th leg of the first tensor has object $(mora[i]), while the second tensor has object $(morb[i])"))
            end
        end
    end

    moradd = Mor(T, mora.objects)

    for key in keys(mora.data)
        moradd[key] = mora[key] + morb[key]
    end

    return moradd
end

function Base.:-(mora::Mor{G, T}, morb::Mor{G, T}) where {T, G<:Group}
    if length(mora.objects)!=length(morb.objects)
        throw(ArgumentError("The first tensor has $(length(mora.objects)) number of legs, while the second tensor has $(length(morb.objects)) number of legs"))
    else
        for i in 1:length(mora.objects)
            if mora[i]!=morb[i]
                throw(ArgumentError("The $i-th leg of the first tensor has object $(mora[i]), while the second tensor has object $(morb[i])"))
            end
        end
    end

    morsubt = Mor(T, mora.objects)

    for key in keys(mora.data)
        morsubt[key] = mora[key] - morb[key]
    end

    return morsubt
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

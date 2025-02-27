# ===========================================
# Abstract Group and Group Element Definitions
# ===========================================
using IterTools

abstract type Group{T} end

"""
CyclicGroup(n)

Construct a structure of group ℤn.

# Example

```
julia> CyclicGroup(4)
CyclicGroup(4)
```
"""
struct CyclicGroup <: Group{Int}
    n::Int
end

"""
DihedralGroup(n)

Construct a structure of group D2n.

# Example

```jldoctest

julia> DihedralGroup(4)
DihedralGroup(4)
```
"""
struct DihedralGroup <: Group{Tuple{Int,Int}}
    n::Int
end

# 辅助函数，用于获取 Group 的值类型
Base.eltype(::Type{<:CyclicGroup}) = Int
Base.eltype(::Type{<:DihedralGroup}) = Tuple{Int,Int}

"""
Definition of the structure GroupElement
"""
struct GroupElement{G<:Group}
    value::eltype(G)  # 自动根据 Group 类型推断 value 的类型
    group::G
end


"""
Construct a group element of group, with value value.

# Parameter:
- value: label of group element
- group: abstract group structure

# Example

```jldoctest
julia> GroupElement(2, CyclicGroup(4))
GroupElement{CyclicGroup}(2, CyclicGroup(4))

julia> GroupElement((2,1), DihedralGroup(4))
GroupElement{DihedralGroup}((2,1), DihedralGroup(4))

```

"""
function GroupElement(value::Any, group::CyclicGroup)
    return GroupElement{CyclicGroup}(mod(value, group.n), group)
end
function GroupElement(value::Tuple, group::DihedralGroup)
    s, r = value
    return GroupElement{DihedralGroup}((mod(s, 2), mod(r, group.n)), group)
end

# ===========================================
# Group Operations
# ===========================================

"""
Generating the tuple of elements of group.

# Example
```
julia> elements(CyclicGroup(3))
(GroupElement{CyclicGroup}(0,CyclicGroup(3)), GroupElement{CyclicGroup}(1,CyclicGroup(3)), GroupElement{CyclicGroup}(2,CyclicGroup(3)))

julia> elements(DihedralGroup(3))
(GroupElement{DihedralGroup}((0,0), DihedralGroup(3)), GroupElement{DihedralGroup}((0,1), DihedralGroup(3)), GroupElement{DihedralGroup}((0,2), DihedralGroup(3)), GroupElement{DihedralGroup}((1,0), DihedralGroup(3)), GroupElement{DihedralGroup}((1,1), DihedralGroup(3)), GroupElement{DihedralGroup}((1,2), DihedralGroup(3)))
```
"""
function elements(group::CyclicGroup)::Tuple
    return ntuple(i -> GroupElement(i - 1, group), group.n)
end
function elements(group::DihedralGroup)::Tuple
    n = group.n
    rotations = ntuple(i -> GroupElement((0, i - 1), group), n)  # (e, r, r², ...)
    reflections = ntuple(i -> GroupElement((1, i - 1), group), n)  # (s, sr, sr², ...)
    return (rotations..., reflections...)  # 将旋转和反射拼接
end

"""
Give the identity element of the group

# Example

```
julia> identity_element(CyclicGroup(3))
GroupElement{CyclicGroup}(0, CyclicGroup(3))
julia> identity_element(DihedralGroup(3))
GroupElement{DihedralGroup}((0,0),DihedralGroup(3))
```

"""
function identity_element(group::CyclicGroup)
    return GroupElement(0, group)
end

function identity_element(group::DihedralGroup)
    return GroupElement((0, 0), group)
end

"""
Give the inverse of the element x

# Example
```
julia> inverse(GroupElement(1, CyclicGroup(3)))
GroupElement{CyclicGroup}(2, CyclicGroup(3))

julia> inverse(GroupElement((0,2), DihedralGroup(3)))
GroupElement{DihedralGroup}((0,1), DihedralGroup(3))
```
"""
function inverse(x::GroupElement{CyclicGroup})
    return GroupElement(-x.value, x.group)
end
function inverse(x::GroupElement{DihedralGroup})
    s, r = x.value
    return GroupElement((-s, (-1)^(s + 1) * r), x.group)
end

"""
Overload the multiplication * to group multiplication.

# Example
```

julia> GroupElement(2, CyclicGroup(3)) * GroupElement(1, CyclicGroup(3))
GroupElement{CyclicGroup}(0, CyclicGroup(3))
```
"""
function Base.:*(x::GroupElement, y::GroupElement)
    if x.group != y.group
        error("Cannot multiply elements from different groups")
    end
    return group_multiply(x, y)
end

# 循环群的乘法
function group_multiply(x::GroupElement{CyclicGroup}, y::GroupElement{CyclicGroup})
    group = x.group
    return GroupElement(x.value + y.value, group)
end

# 二面体群的乘法
function group_multiply(x::GroupElement{DihedralGroup}, y::GroupElement{DihedralGroup})
    group = x.group
    s1, r1 = x.value
    s2, r2 = y.value
    return GroupElement((mod(s1 + s2, 2), mod((-1)^s2 * r1 + r2, group.n)), group)
end

"""
Multiplying a tuple of group elements.

# Example

```
julia> multiply((GroupElement(1, CyclicGroup(4)), GroupElement(2, CyclicGroup(4))))
GroupElement{CyclicGroup}(3, CyclicGroup(4))
```
"""
function multiply(elements::Tuple{Vararg{GroupElement}})
    if length(elements) == 0
        error("Cannot multiply an empty tuple of group elements")
    elseif length(elements) == 1
        return elements[1]
    else
        result = elements[1]
        for i in 2:length(elements)
            result *= elements[i]
        end
        return result
    end
end

"""
Overload the == to group case.
"""
function Base.:(==)(x::GroupElement{CyclicGroup}, y::GroupElement{CyclicGroup})
    return x.value == y.value
end

function Base.:(==)(x::GroupElement{DihedralGroup}, y::GroupElement{DihedralGroup})
    return x.value == y.value
end

"""
```group_tree(g, n)```

Generating an iterator of n-tuple of group elements, such that all elements are multiplied to g.

# Input Parameter:

- g: a group element
- n: the number of group elements in a tuple

# Example

```
julia> x = group_tree(GroupElement(0, CyclicGroup(2)), 2)

Base.Generator{Base.Iterators.ProductIterator{Tuple{Tuple{GroupElement{CyclicGroup}, GroupElement{CyclicGroup}}}}, var"#18#20"{GroupElement{CyclicGroup}}}(var"#18#20"{GroupElement{CyclicGroup}}(GroupElement{CyclicGroup}(0, CyclicGroup(2))), Base.Iterators.ProductIterator{Tuple{Tuple{GroupElement{CyclicGroup}, GroupElement{CyclicGroup}}}}(((GroupElement{CyclicGroup}(0, CyclicGroup(2)), GroupElement{CyclicGroup}(1, CyclicGroup(2))),)))

julia> for i in x
            println(i)
       end

(GroupElement{CyclicGroup}(0, CyclicGroup(2)), GroupElement{CyclicGroup}(0, CyclicGroup(2)))
(GroupElement{CyclicGroup}(1, CyclicGroup(2)), GroupElement{CyclicGroup}(1, CyclicGroup(2)))
```
"""
function group_tree(g::GroupElement, n::Int)
    if n < 1
        error("Invalid length of the tuple")
    elseif n == 1
        return ((g,),)
    else
        group = g.group  # 从 g 中推断所属群
        elem = elements(group)
        tup = IterTools.product((elem for _ in 1:(n-1))...)
        return Iterators.map(x -> (x..., inverse(inverse(g) * multiply(x))), tup)
    end
end

"""
Generating a iterator, whose terms are n-tuple of group elements.

# Example

```
julia>  for i in group_iter(CyclicGroup(2), 2)
            println(i)
        end
(GroupElement{CyclicGroup}(0, CyclicGroup(2)), GroupElement{CyclicGroup}(0, CyclicGroup(2)))
(GroupElement{CyclicGroup}(1, CyclicGroup(2)), GroupElement{CyclicGroup}(0, CyclicGroup(2)))
(GroupElement{CyclicGroup}(0, CyclicGroup(2)), GroupElement{CyclicGroup}(1, CyclicGroup(2)))
(GroupElement{CyclicGroup}(1, CyclicGroup(2)), GroupElement{CyclicGroup}(1, CyclicGroup(2)))
```
"""
function group_iter(group::Group, n::Int)
    elem = elements(group)
    return IterTools.product((elem for _ in 1:n)...)
end

"""
Helper to verify group axioms.

# Input:
- g: group

"""
function verify_group_axioms(g::Group)
    elts = elements(g)
    println("Verifying Group Axioms for ", typeof(g))

    # Identity element
    id = identity_element(g)
    println("Identity: ", all(x -> (x * id).value == x.value && (id * x).value == x.value, elts))

    # Associativity
    println("Associativity: ", all(a -> all(b -> all(c -> ((a * b) * c).value == (a * (b * c)).value, elts), elts), elts))

    # Inverses
    println("Inverses: ", all(x -> (x * inverse(x)).value == id.value && (inverse(x) * x).value == id.value, elts))
end


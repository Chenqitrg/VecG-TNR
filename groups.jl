# ===========================================
# Abstract Group and Group Element Definitions
# ===========================================
using IterTools

abstract type Group{T} end

struct CyclicGroup <: Group{Int}
    n::Int
end

struct DihedralGroup <: Group{Tuple{Int, Int}}
    n::Int
end

# 辅助函数，用于获取 Group 的值类型
Base.eltype(::Type{<:CyclicGroup}) = Int
Base.eltype(::Type{<:DihedralGroup}) = Tuple{Int, Int}

# 群元素定义
struct GroupElement{G<:Group}
    value::eltype(G)  # 自动根据 Group 类型推断 value 的类型
    group::G
end


# 群元素构造函数，支持循环群和二面体群
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

# 群元素生成函数，多重分派优化
function elements(group::CyclicGroup)::Tuple
    return ntuple(i -> GroupElement(i - 1, group), group.n)
end

function elements(group::DihedralGroup)::Tuple
    n = group.n
    rotations = ntuple(i -> GroupElement((0, i - 1), group), n)  # (e, r, r², ...)
    reflections = ntuple(i -> GroupElement((1, i - 1), group), n)  # (s, sr, sr², ...)
    return (rotations..., reflections...)  # 将旋转和反射拼接
end

# 单位元
function identity(group::CyclicGroup)
    return GroupElement(0, group)
end

function identity(group::DihedralGroup)
    return GroupElement((0, 0), group)
end

# 逆元
function inverse(x::GroupElement{CyclicGroup})
    return GroupElement(-x.value, x.group)
end

function inverse(x::GroupElement{DihedralGroup})
    s, r = x.value
    return GroupElement((-s, (-1)^(s + 1) * r), x.group)
end

# 通用乘法接口
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

function Base.:(==)(x::GroupElement{CyclicGroup}, y::GroupElement{CyclicGroup})
    return x.value == y.value
end

function Base.:(==)(x::GroupElement{DihedralGroup}, y::GroupElement{DihedralGroup})
    return x.value == y.value
end

function group_tree(g::GroupElement, n::Int)
    if n < 1
        error("Invalid length of the tuple")
    elseif n == 1
        return ((g,),)
    else
        group = g.group  # 从 g 中推断所属群
        elem = elements(group)
        tup = IterTools.product((elem for _ in 1:(n - 1))...)
        return Iterators.map(x -> (x..., inverse(inverse(g) * multiply(x))), tup)
    end
end

# Helper to verify group axioms
function verify_group_axioms(g::Group)
    elts = elements(g)
    println("Verifying Group Axioms for ", typeof(g))

    # Identity element
    id = identity(g)
    println("Identity: ", all(x -> (x * id).value == x.value && (id * x).value == x.value, elts))

    # Associativity
    println("Associativity: ", all(a -> all(b -> all(c -> ((a * b) * c).value == (a * (b * c)).value, elts), elts), elts))

    # Inverses
    println("Inverses: ", all(x -> (x * inverse(x)).value == id.value && (inverse(x) * x).value == id.value, elts))
end


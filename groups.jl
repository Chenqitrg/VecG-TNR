# ===========================================
# Abstract Group and Group Element Definitions
# ===========================================
using IterTools

abstract type Group end

struct CyclicGroup <: Group
    n::Int
end

struct DihedralGroup <: Group
    n::Int
end

struct GroupElement{G<:Group, T}
    value::T
    group::G
end

# 群元素构造函数，支持循环群和二面体群
function GroupElement(value::Any, group::G) where {G<:Group}
    if group isa CyclicGroup
        return GroupElement{G, Int}(mod(value, group.n), group)
    elseif group isa DihedralGroup
        s, r = value
        return GroupElement{G, Tuple{Int, Int}}((mod(s, 2), mod(r, group.n)), group)
    end
    error("Unsupported group type")
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
function inverse(x::GroupElement{CyclicGroup, Int})
    return GroupElement(-x.value, x.group)
end

function inverse(x::GroupElement{DihedralGroup, Tuple{Int, Int}})
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
function group_multiply(x::GroupElement{CyclicGroup, Int}, y::GroupElement{CyclicGroup, Int})
    group = x.group
    return GroupElement(x.value + y.value, group)
end

# 二面体群的乘法
function group_multiply(x::GroupElement{DihedralGroup, Tuple{Int, Int}}, y::GroupElement{DihedralGroup, Tuple{Int, Int}})
    group = x.group
    s1, r1 = x.value
    s2, r2 = y.value
    return GroupElement((mod(s1 + s2, 2), mod((-1)^s2 * r1 + r2, group.n)), group)
end

function multiply(elements::Tuple{Vararg{GroupElement}})
    if length(elements) == 0
        error("Cannot multiply an empty tuple of group elements")
    end
    result = elements[1]
    for i in 2:length(elements)
        result *= elements[i]
    end
    return result
end

function group_tree(group::G, g::GroupElement{G}, n::Int) where G<:Group
    if n < 1
        error("Invalid length of the tuple")
    elseif n == 1
        return (g,)
    else
        elem = elements(group)
        tup = IterTools.product((elem for _ in 1:(n-1))...)
        return Iterators.map(x -> (x..., inverse(inverse(g) * multiply(x))), tup)
    end
end

# ===========================================
# Display Functions
# ===========================================

# 群元素的显示函数
function Base.show(io::IO, x::GroupElement{CyclicGroup, Int})
    if x.value == 0
        print(io, "e")
    elseif x.value == 1
        print(io, "a")
    else
        print(io, "a", superscript(x.value))
    end
end

function Base.show(io::IO, x::GroupElement{DihedralGroup, Tuple{Int, Int}})
    s, r = x.value
    if s == 0 && r == 0
        print(io, "e")
    elseif s == 0 && r == 1
        print(io, "r")
    elseif s == 0 && r > 0
        print(io, "r", superscript(r))
    elseif s == 1 && r == 0
        print(io, "s")
    elseif s == 1 && r == 1
        print(io, "sr")
    else
        print(io, "sr", superscript(r))
    end
end

function Base.show(io::IO, G::Group)
    if G isa CyclicGroup
        print(io, "ℤ", subscript(G.n))
    elseif G isa DihedralGroup
        print(io, "𝔻", subscript(2 * G.n))
    end
end

function subscript(n::Integer)
    # 定义一个映射，将字符 '0'-'9' 转换为 Unicode 下角标字符
    subs = Dict('0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄',
                '5' => '₅', '6' => '₆', '7' => '₇', '8' => '₈', '9' => '₉')
    # 将数字 n 转换为字符串，逐位映射到下角标字符，并连接成新的字符串
    return join([subs[c] for c in string(n)])
end

function superscript(n::Integer)
    # 定义一个映射，将字符 '0'-'9' 转换为 Unicode 下角标字符
    subs = Dict('0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
                '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹')
    # 将数字 n 转换为字符串，逐位映射到下角标字符，并连接成新的字符串
    return join([subs[c] for c in string(n)])
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


# ===========================================
# Abstract Group and Group Element Definitions
# ===========================================

abstract type Group end

struct CyclicGroup <: Group
    n::Int
end

struct DihedralGroup <: Group
    n::Int
end

struct GroupElement{G<:Group}
    value::Any
    group::G
    function GroupElement(value::Any, group::G) where {G<:Group}
        if group isa CyclicGroup  # Check the instance, not the type parameter
            return new{G}(mod(value, group.n), group)
        elseif group isa DihedralGroup
            s, r = value
            return new{G}((mod(s,2), mod(r,group.n)), group)
        end
        error("Unsupported group type")
    end
end

function elements(group::Group)
    if group isa CyclicGroup  # Check the instance, not the type parameter
        return [GroupElement(i, group) for i in 0:(group.n - 1)]
    elseif group isa DihedralGroup
        return vec([GroupElement((s, r), group) for r in 0:(group.n - 1), s in 0:1])
    end
    error("Unsupported group type")
    return 
end

function Base.show(io::IO, x::GroupElement)
    if x isa GroupElement{DihedralGroup}
        s, r = x.value
        if s == 0 && r == 0
            print(io, "e")  # 单位元显示为 e
        elseif s == 0
            print(io, "r^", r)  # 旋转元素显示为 r^k
        elseif s == 1 && r == 0
            print(io, "s")  # 反射元素显示为 s
        else
            print(io, "sr^", r)  # 反射+旋转显示为 sr^k
        end
    elseif x isa GroupElement{CyclicGroup}
        if x.value == 0
            print(io, "e")
        else
            print(io, "a^", x.value)
        end
    end
end

# @show elements(CyclicGroup(5))

function identity(group::Group)
    if group isa CyclicGroup
        return GroupElement(0, group)
    elseif group isa DihedralGroup
        return GroupElement((0,0),group)
    end
    error("Unsupported group type")
    return 
end

identity(DihedralGroup(3))

function inverse(x::GroupElement)
    if x isa GroupElement{CyclicGroup}
        return GroupElement(-x.value, x.group)
    elseif x isa GroupElement{DihedralGroup}
        s, r = x.value
        return GroupElement((-s, (-1)^(s+1) * r), x.group)
    end
    error("Unsupported group type")
    return 
end

# inverse(GroupElement((1,0), DihedralGroup(3)))

function Base.:*(x::GroupElement, y::GroupElement)
    if x.group == y.group
        group = x.group
        if x isa GroupElement{CyclicGroup}
            return GroupElement(x.value + y.value, group)
        elseif x isa GroupElement{DihedralGroup}
            s1, r1 = x.value
            s2, r2 = y.value
            return GroupElement((mod(s1 + s2, 2), mod((-1)^s2 * r1 + r2, group.n)), group)
        end
        error("Unsupported group type")
    end
    error("Not allowed to multiply")
    return 
end

function multiply(v::Vector{GroupElement{G}}) where G <: Group
    res = v[1]
    for j = 2 : length(v)
        res = res * v[j]
    end
    return res
end
# x = GroupElement((3,1), DihedralGroup(4))
# y = GroupElement((1,2), DihedralGroup(4))
# w = GroupElement((0,1), DihedralGroup(4))
# v = [x, y, w]
# @show multiply(v)

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

# verify_group_axioms(DihedralGroup(5))


# # ===========================================
# # Direct Product Group G x H
# # ===========================================

# struct DirectProductGroup{G<:Group, H<:Group} <: Group
#     g::G
#     h::H
# end

# function identity(g::DirectProductGroup)
#     return GroupElement((identity(g.g), identity(g.h)), g)
# end

# function inv(x::GroupElement{DirectProductGroup})
#     g = x.group
#     gx, hx = x.value
#     return GroupElement((inv(gx), inv(hx)), g)
# end

# function *(x::GroupElement{DirectProductGroup}, y::GroupElement{DirectProductGroup})
#     g = x.group
#     gx1, hx1 = x.value
#     gx2, hx2 = y.value
#     return GroupElement((gx1 * gx2, hx1 * hx2), g)
# end



# function elements(g::DirectProductGroup)
#     return [GroupElement((gx, hy), g) for gx in elements(g.g), hy in elements(g.h)]
# end

# # ===========================================
# # Example Usage and Testing
# # ===========================================

# # Test Zn (e.g., Z4 and Z5)
# Z4 = CyclicGroup(4)
# println("Elements of Z4: ", elements(Z4))
# verify_group_axioms(Z4)

# # Test D2n (e.g., D8 for n=4)
# D8 = DihedralGroup(4)
# println("\nElements of D8: ", elements(D8))
# verify_group_axioms(D8)

# # Test Direct Product G x H (e.g., Z2 x Z3)
# Z2 = CyclicGroup(2)
# Z3 = CyclicGroup(3)
# Z2xZ3 = DirectProductGroup(Z2, Z3)
# println("\nElements of Z2 x Z3: ", elements(Z2xZ3))
# verify_group_axioms(Z2xZ3)
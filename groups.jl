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
            return new{G}((mod(s,2), mod(r,4)), group)
        end
        error("Unsupported group type")
    end
end

D4 = DihedralGroup(4)
el = GroupElement((3,4), D4)
println(el)

# function elements(g::CyclicGroup)
#     return [GroupElement(i, g) for i in 0:(g.n - 1)]
# end



# # Group operations
# identity(g::Group) = throw(NotImplementedError("Identity not defined for this group"))
# inv(x::GroupElement) = throw(NotImplementedError("Inverse not defined for this group"))
# *(x::GroupElement, y::GroupElement) = throw(NotImplementedError("Operation not defined for this group"))
# elements(g::Group) = throw(NotImplementedError("Elements not defined for this group"))

# # Helper to verify group axioms
# function verify_group_axioms(g::Group)
#     elts = elements(g)
#     println("Verifying Group Axioms for ", typeof(g))

#     # Identity element
#     id = identity(g)
#     println("Identity: ", all(x -> (x * id).value == x.value && (id * x).value == x.value, elts))

#     # Associativity
#     println("Associativity: ", all(a -> all(b -> all(c -> ((a * b) * c).value == (a * (b * c)).value, elts), elts), elts))

#     # Inverses
#     println("Inverses: ", all(x -> (x * inv(x)).value == id.value && (inv(x) * x).value == id.value, elts))
# end

# ===========================================
# Cyclic Group Zn
# ===========================================

# struct CyclicGroup <: Group
#     n::Int
# end


# function identity(g::CyclicGroup)
#     return GroupElement(0, g)
# end

# function inv(x::GroupElement{CyclicGroup})
#     g = x.group
#     return GroupElement(mod(-x.value, g.n), g)
# end

# function Base.:*(x::GroupElement{CyclicGroup}, y::GroupElement{CyclicGroup})
#     g = x.group
#     @assert g === y.group "Elements must belong to the same group"
#     return GroupElement(mod(x.value + y.value, g.n), g)
# end

# function elements(g::CyclicGroup)
#     return [GroupElement(i, g) for i in 0:(g.n - 1)]
# end

# Z3 = CyclicGroup(3)
# elements(Z3)
# ===========================================
# Dihedral Group D2n
# ===========================================



# function identity(g::DihedralGroup)
#     return GroupElement((0, 0), g)  # (reflection, rotation)
# end

# function inv(x::GroupElement{DihedralGroup})
#     g = x.group
#     s, r = x.value
#     if s == 0  # rotation
#         return GroupElement((0, mod(-r, g.n)), g)
#     else  # reflection
#         return GroupElement((1, mod(r, g.n)), g)
#     end
# end

# function Base.:*(x::GroupElement{DihedralGroup}, y::GroupElement{DihedralGroup})
#     g = x.group
#     @assert g === y.group "Elements must belong to the same group"
#     (s1, r1) = x.value
#     (s2, r2) = y.value
#     return GroupElement((mod(s1 + s2, 2), mod((-1)^s2 * r1 + r2, g.n)), g)
# end

# function elements(g::DihedralGroup)
#     return vec([GroupElement((s, r), g) for r in 0:(g.n - 1), s in 0:1])
# end

# function Base.show(io::IO, x::GroupElement{DihedralGroup})
#     s, r = x.value
#     if s == 0 && r == 0
#         print(io, "e")  # 单位元显示为 e
#     elseif s == 0
#         print(io, "r^", r)  # 旋转元素显示为 r^k
#     elseif s == 1 && r == 0
#         print(io, "s")  # 反射元素显示为 s
#     else
#         print(io, "sr^", r)  # 反射+旋转显示为 sr^k
#     end
# end

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
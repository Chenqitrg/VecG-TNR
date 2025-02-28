using TensorKit

struct VecD8 <: Sector
    s::Int
    r::Int
end

Base.one(::Type{VecD8}) = VecD8(0,0)
Base.conj(g::VecD8) = VecD8(mod(-g.s,2), mod((-1)^(g.s + 1) * g.r,4))
dim(g::VecD8) = 1
Base.isreal(::Type{VecD8}) = true
TensorKit.FusionStyle(::Type{VecD8}) = UniqueFusion()
TensorKit.BraidingStyle(::Type{VecD8}) = Bosonic()
TensorKit.Nsymbol(a::VecD8, b::VecD8, c::VecD8) = (mod(a.s + b.s, 2), mod((-1)^b.s * a.r + b.r, 4)) == (c.s, c.r) ? 1 : 0
⊗(a::VecD8, b::VecD8) = VecD8(mod(a.s + b.s, 2), mod((-1)^b.s * a.r + b.r, 4))
Base.isless(a::VecD8, b::VecD8) = a.s < b.s || (a.s == b.s && a.r < b.r)
function Base.iterate(::TensorKit.SectorValues{VecD8},s=0,r=0)
    if s==0 && r < 4
        return (VecD8(s,r), (s, r+1))
    elseif s==0 && r==4
        return (VecD8(0,4), (1,0))
    elseif s==1 && r < 4
        return (VecD8(s,r), (s, r+1))
    else
        return nothing
    end
end
TensorKit.Fsymbol(a::VecD8, b::VecD8, c::VecD8, d::VecD8, e::VecD8, f::VecD8) = 1
Base.length(::TensorKit.SectorValues{VecD8}) = 8
Base.getindex(::TensorKit.SectorValues{VecD8}, i::Int) = VecD8((i-1)÷4, (i-1)%4)
TensorKit.findindex(::TensorKit.SectorValues{VecD8}, c::VecD8) = c.r + 4 * c.s + 1
V1 = Vect[VecD8](VecD8(0,0) => 1, VecD8(1,2) => 2)
V2 = Vect[VecD8](VecD8(0,0) => 1, VecD8(1,2) => 2)
V3 = Vect[VecD8](VecD8(0,0) => 1, VecD8(1,2) => 2)
V4 = Vect[VecD8](VecD8(0,0) => 1, VecD8(1,2) => 2)

function ⊗(V1::GradedSpace{VecD8, NTuple{8, Int64}}, V2::GradedSpace{VecD8, NTuple{8, Int64}})
    D8 = TensorKit.SectorValues{VecD8}()
    multiplicity = zeros(Int, 8)
    for i in 1:8, j in 1:8
        delta = V1.dims[i] * V2.dims[j]
        g = D8[i] ⊗ D8[j]
        multiplicity[TensorKit.findindex(D8, g)] += delta
    end
    return GradedSpace{VecD8, NTuple{8, Int64}}(Tuple(multiplicity), false)
end

A = rand(V1 ⊗ V4 ← V2 ⊗ V3)
VecD8(1,2) ⊗ VecD8(0,0)
GradedSpace{VecD8, NTuple{8, Int64}}((1, 0, 0, 0, 0, 0, 0, 0), false)

# for a in sectors(V1), b in sectors(V2)
#     @show a ⊗ b
#     # for c in a ⊗ b
#     #     println(c)
#     # end
# end
# V1 ⊗ V2
# function fuse(V₁::GradedSpace{I}, V₂::GradedSpace{I}) where {I<:Sector}
#     dims = SectorDict{I,Int}()
#     for a in sectors(V₁), b in sectors(V₂)
#         for c in a ⊗ b
#             dims[c] = get(dims, c, 0) + Nsymbol(a, b, c) * dim(V₁, a) * dim(V₂, b)
#         end
#     end
#     return typeof(V₁)(dims)
# end
# fuse(V1, V2)

# sectors = TensorKit.SectorValues{VecD8}()
# println(sectors[1])  # 输出 VecD8(0, 0)
# println(sectors[5])  # 输出 VecD8(1, 0)


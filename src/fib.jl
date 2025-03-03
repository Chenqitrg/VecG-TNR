# abstract type UFC end
# struct OneDFusionSpaceUFC <: UFC
#     order::Int
#     fusion_rule::Dict{Tuple{Int, Int, Int}, Int}
#     F::Dict{Tuple{Int, Int, Int, Int}, Array{Float64}}
# end

# phi = (sqrt(5)+1)/2

# fib_fusion_rule = Dict(
#     (0, 0, 0) => 1,
#     (0, 1, 1) => 1,
#     (1, 0, 1) => 1,
#     (1, 1, 0) => 1,
#     (1, 1, 1) => 1
# )

# # 定义 Fibonacci category 的六角系数 (F-symbols)
# fib_F = Dict(
#     (0,0,0,0)=>[1.],
#     (0,0,1,1)=>[1.],
#     (0,1,0,1)=>[1.],
#     (0,1,1,0)=>[1.],
#     (0,1,1,1)=>[1.],
#     (1,0,0,1)=>[1.],
#     (1,0,1,0)=>[1.],
#     (1,0,1,1)=>[1.],
#     (1,1,0,0)=>[1.],
#     (1,1,0,1)=>[1.],
#     (1,1,1,1)=>[1/phi 1/sqrt(phi); 1/sqrt(phi) -1/phi]
# )

# Fib = OneDFusionSpaceUFC(2, fib_fusion_rule, fib_F)

# struct Irr
#     label::Int
#     cat::OneDFusionSpaceUFC
#     function Irr(label::Int, cat::OneDFusionSpaceUFC)
#         n = cat.order
#         return new(mod(label, n), cat)
#     end
# end

# function get_order(x::Irr)
#     cat = x.cat
#     order = cat.order
#     return order
# end

# function sobs(cat::OneDFusionSpaceUFC)
#     order = cat.order
#     tup = ()
#     for i in 0:order-1
#         tup = (tup..., Irr(i, cat))
#     end
#     return tup
# end

# struct Obj
#     sumd::Dict{Irr, Int}
# end

# function Obj(pairs::Pair{Irr, Int}...)
#     order = get_order(pairs[1].first)
#     cat = pairs[1].first.cat
#     sumd = Dict{Irr, Int}()
#     for pair in pairs
#         sumd[pair.first] = pair.second
#     end

#     for i in sobs(cat)
#         if !haskey(sumd, i)
#             sumd[i] = 0
#         end
#     end

#     return Obj(sumd)
# end

# function get_order(obj::Obj)
#     return length(obj.sumd)
# end

# function get_cat(obj::Obj)
#     k = keys(obj.sumd)
#     key1 = first(k)
#     return key1.cat
# end

# function Base.show(io::IO, x::Irr)
#     if x.cat == Fib
#         if x.label == 0
#             print(io, "𝕀")
#         elseif x.label == 1
#             print(io, "τ")
#         end
#     end
# end

# # Custom display for Object
# function Base.show(io::IO, x::Obj)
#     sumd = x.sumd
    
#     if sumd[i] != 0
#         if op == true
#             print(io, " ⊕ ")
#         end
#         if sumd[i] == 1
#             print(io, i)
#         else
#             print(io, sumd[i], i)
#         end
#         op = true
#     end
# end

e = Irr(0, Fib)
t = Irr(1, Fib)
sobs(Fib)

ob = Obj(e=>2)

get_cat(ob)
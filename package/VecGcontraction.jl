# include("groups.jl")
# include("VecGtensor.jl")
# include("display.jl")
using LinearAlgebra
function VecG_tensordot(A::Mor{G,T}, B::Mor{G,T}, leg_cont::Int) where {T, G<:Group}
    A_legs = length(A.objects)
    B_legs = length(B.objects)
    A_legs_remain = A_legs-leg_cont
    B_legs_remain = B_legs - leg_cont

    group = get_group(A)
    for i in 1:leg_cont
        if A[end-i+1]!=dual_obj(B[i])
            throw(ArgumentError("Contracted object $(A[end-i+1]) is not the dual of $(B[i])"))
        end
    end

    newobjs = (A[1:A_legs_remain]..., B[leg_cont+1:B_legs]...)

    Cont = Mor(T, newobjs)

    for g_bridge in elements(group)
        for sect_A in group_tree(g_bridge, A_legs_remain), sect_B in group_tree(inverse(g_bridge), B_legs_remain)
            tot_sect = (sect_A...,sect_B...)
            tot_sect_size = get_sector_size(Cont, tot_sect)
            Cont[tot_sect...] = zeros(T,tot_sect_size)
            for sect_cont in group_tree(g_bridge, leg_cont)
                sect_cont_prime = inverse.(reverse(sect_cont))
                A_tot_sect = (sect_A..., sect_cont_prime...)
                B_tot_sect = (sect_cont..., sect_B...)
                A_block = A[A_tot_sect...]
                B_block = B[B_tot_sect...]
                A_multiplicity = get_sector_size(A, A_tot_sect)
                B_multiplicity = get_sector_size(B, B_tot_sect)
                A_dim_left = prod(A_multiplicity[1:A_legs_remain])
                A_dim_cont = prod(A_multiplicity[A_legs_remain+1:A_legs])
                B_dim_cont = prod(B_multiplicity[1:leg_cont])
                B_dim_left = prod(B_multiplicity[leg_cont+1:B_legs])
                Temp_mat = reshape(reshape(A_block, (A_dim_left, A_dim_cont)) * reshape(B_block, (B_dim_cont,B_dim_left)), tot_sect_size)
                Cont[tot_sect...] = Cont[tot_sect...] + Temp_mat
            end
        end
    end

    return Cont
end

function VecG_tensordot(A::Mor{G,T}, B::Mor{G,T}, A_cont::Tuple{Vararg{Int}}, B_cont::Tuple{Vararg{Int}}) where {T, G<:Group}
    dimA = length(A.objects)
    dimB = length(B.objects)
    if !(is_accend(A_cont, dimA)&&is_accend(B_cont, dimB)) || length(A_cont)!=length(B_cont)
        throw(ArgumentError("The contraction indices $A_cont and $B_cont are illegal"))
    end
    leg_cont = length(A_cont)

    A_perm = (A_cont[end]-dimA+1):A_cont[end]
    A_perm = Tuple(map(x->mod(x-1,dimA)+1, A_perm))
    B_perm = to_perm(B_cont, dimB)

    newA = VecG_permutedims(A, A_perm)
    newB = VecG_permutedims(B, B_perm)

    return VecG_tensordot(newA, newB, leg_cont)
end


function VecG_factorize(mor::Mor, n_leg_split::Tuple{Vararg{Int}}, Dcut::Int, epsilon::Float64)
    U, S, V = VecG_svd(mor, n_leg_split, Dcut, epsilon)
    group = get_group(S)
    for k in elements(group)
        S[k, inverse(k)] = sqrt.(S[k, inverse(k)])
    end
    F = VecG_tensordot(U, S, 1)
    K = VecG_tensordot(V, VecG_permutedims(S, (2,1)), 1)
    
    return F, K
end

function VecG_partial_trace(mor::Mor{G, T}, leg_cont::Int) where {T, G<:Group}
    group = get_group(mor)
    e = identity(group)
    n_leg = length(mor.objects)

    if 2*leg_cont > n_leg
        throw(ArgumentError("The number of traced legs $(2*leg_cont) is larger than the number of total legs $n_leg"))
    end

    for i = 1 : leg_cont
        if mor[n_leg-2*leg_cont+i] != dual_obj(mor[n_leg-i+1])
            throw(ArgumentError("The $(n_leg-2*leg_cont+i)-th leg is $(mor[n_leg-2*leg_cont+i]), not the dual of $(n_leg-i+1)-th leg $(mor[n_leg-i+1])"))
        end
    end

    leg_remain = n_leg - 2 * leg_cont

    objs = mor[1:leg_remain]

    new_mor = zero_mor(T, objs)
    for sect_remain in group_tree(e, leg_remain)
        for sect_tr in group_iter(group, leg_cont)
            inv_sect_tr = reverse(inverse.(sect_tr))
            new_mor[sect_remain...] = new_mor[sect_remain...] + partial_trace(mor[sect_remain..., sect_tr..., inv_sect_tr...], leg_cont)
        end
    end

    return new_mor
end

# Z4 = CyclicGroup(4)
# e = GroupElement(0, Z4)
# a = GroupElement(1, Z4)

# A = Obj(e=>2, a=>2, a*a=>3, a*a*a=>4)
# B = Obj(e=>2, a=>4, a*a=>3, a*a*a=>2)

# A_dual = dual_obj(A)
# B_dual = dual_obj(B)

# T1 = random_mor(Float64, (A,B,A,B))
# T2 = random_mor(Float64, (B_dual, A_dual, B_dual,B))

# # T = VecG_tensordot(T1, T2, 2)
# Tp = VecG_tensordot(T1, T2, (2,3), (2,3))

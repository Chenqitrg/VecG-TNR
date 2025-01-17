
function VecG_tensordot(A::Mor{G,T}, B::Mor{G,T}, leg_cont::Int) where {T,G<:Group}
    A_legs = length(A.objects)
    B_legs = length(B.objects)
    A_legs_remain = A_legs - leg_cont
    B_legs_remain = B_legs - leg_cont

    group = get_group(A)
    for i in 1:leg_cont
        if A[end-i+1] != dual_obj(B[i])
            throw(ArgumentError("Contracted object $(A[end-i+1]) is not the dual of $(B[i])"))
        end
    end

    if (A_legs_remain + B_legs_remain) == 0
        Cont = 0.
        for sect_cont in group_tree(identity_element(group), leg_cont)
            sect_cont_prime = inverse.(reverse(sect_cont))
            A_block = A[sect_cont_prime...]
            B_block = B[sect_cont...]
            A_multiplicity = get_sector_size(A, sect_cont_prime)
            B_multiplicity = get_sector_size(B, sect_cont)
            A_dim_cont = prod(A_multiplicity)
            B_dim_cont = prod(B_multiplicity)
            A_reshaped = reshape(A_block, A_dim_cont)
            B_reshaped = reshape(permutedims(B_block, reverse(1:leg_cont)), B_dim_cont) # Be careful of the reshaping order: the reshaping needs to be compatible with the contraction order
            Temp_mat = transpose(A_reshaped) * B_reshaped
            Cont = Cont + Temp_mat[1]
        end
        return Cont
    else
        newobjs = (A[1:A_legs_remain]..., B[leg_cont+1:B_legs]...)

        Cont = Mor(T, newobjs)

        for g_bridge in elements(group)
            # @show g_bridge
            for sect_A in group_tree(g_bridge, A_legs_remain), sect_B in group_tree(inverse(g_bridge), B_legs_remain)
                # @show sect_A
                # @show sect_B
                tot_sect = (sect_A..., sect_B...)
                tot_sect_size = get_sector_size(Cont, tot_sect)
                Cont[tot_sect...] = zeros(T, tot_sect_size)
                for sect_cont in group_tree(g_bridge, leg_cont)
                    sect_cont_prime = inverse.(reverse(sect_cont))
                    A_tot_sect = (sect_A..., sect_cont_prime...)
                    B_tot_sect = (sect_cont..., sect_B...)
                    # @show A_tot_sect
                    # @show B_tot_sect
                    A_block = A[A_tot_sect...]
                    B_block = B[B_tot_sect...]
                    A_multiplicity = get_sector_size(A, A_tot_sect)
                    B_multiplicity = get_sector_size(B, B_tot_sect)
                    A_dim_left = prod(A_multiplicity[1:A_legs_remain])
                    A_dim_cont = prod(A_multiplicity[A_legs_remain+1:A_legs])
                    B_dim_cont = prod(B_multiplicity[1:leg_cont])
                    B_dim_left = prod(B_multiplicity[leg_cont+1:B_legs])

                    A_reshaped = reshape(A_block, (A_dim_left, A_dim_cont))
                    B_reshaped = reshape(permutedims(B_block, (reverse(1:leg_cont)..., leg_cont+1:B_legs...)), (B_dim_cont, B_dim_left)) # Be careful of the reshaping order: the reshaping needs to be compatible with the contraction order
                    Temp_mat = reshape(A_reshaped * B_reshaped, tot_sect_size)
                    Cont[tot_sect...] = Cont[tot_sect...] + Temp_mat
                end
            end
        end

        return Cont
    end
end

function VecG_tensordot(A::Mor{G,T}, B::Mor{G,T}, A_cont::Tuple{Vararg{Int}}, B_cont::Tuple{Vararg{Int}}) where {T,G<:Group}
    dimA = length(A.objects)
    dimB = length(B.objects)
    if !(is_accend(A_cont, dimA) && is_decend(B_cont, dimB)) || length(A_cont) != length(B_cont)
        throw(ArgumentError("The contraction indices $A_cont and $B_cont are illegal"))
    end
    leg_cont = length(A_cont)

    A_perm = (A_cont[end]-dimA+1):A_cont[end]
    B_perm = B_cont[end]:(B_cont[end]+dimB-1)
    A_perm = Tuple(map(x -> mod(x - 1, dimA) + 1, A_perm))
    B_perm = Tuple(map(x -> mod(x - 1, dimB) + 1, B_perm))

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
    K = VecG_tensordot(V, VecG_permutedims(S, (2, 1)), 1)

    return F, K
end

function VecG_partial_trace(mor::Mor{G,T}, leg_cont::Int) where {T,G<:Group}
    group = get_group(mor)
    e = identity_element(group)
    n_leg = length(mor.objects)

    if 2 * leg_cont > n_leg
        throw(ArgumentError("The number of traced legs $(2*leg_cont) is larger than the number of total legs $n_leg"))
    end

    for i = 1:leg_cont
        if mor[n_leg-2*leg_cont+i] != dual_obj(mor[n_leg-i+1])
            throw(ArgumentError("The $(n_leg-2*leg_cont+i)-th leg is $(mor[n_leg-2*leg_cont+i]), not the dual of $(n_leg-i+1)-th leg $(mor[n_leg-i+1])"))
        end
    end

    if 2 * leg_cont == n_leg
        num = 0.0
        for sect_tr in group_iter(group, leg_cont)
            inv_sect_tr = reverse(inverse.(sect_tr))
            num = num + reshape(partial_trace(mor[sect_tr..., inv_sect_tr...], leg_cont), 1)[1]
        end
        return num
    else
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
end

function VecG_square(mor::Mor{G,T}) where {T,G<:Group}
    n = length(mor.objects)
    return VecG_tensordot(mor, VecG_dag(mor), n)
end

function VecG_outprod(mor1::Mor{G,T}, mor2::Mor{G,T}) where {T,G<:Group}
    group = get_group(mor1)
    legs1 = length(mor1.objects)
    legs2 = length(mor2.objects)
    e = identity_element(group)
    objs = (mor1.objects..., mor2.objects...)
    newmor = Mor(T, objs)
    for sect in group_tree(e, legs1+legs2)
        sizeT = get_sector_size(newmor, sect)
        if multiply(sect[1:legs1])!=identity_element(group) || multiply(sect[legs1+1:legs1+legs2])!=identity_element(group)
            newmor[sect...] = zeros(T, sizeT)
        else
            sect1 = sect[1:legs1]
            sect2 = sect[legs1+1:legs1+legs2]
            temp = outer_product(mor1[sect1...], mor2[sect2...])
            newmor[sect...] = temp
        end
    end
    return newmor
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

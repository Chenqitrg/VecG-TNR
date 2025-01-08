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
            Cont[tot_sect...] = zeros(tot_sect_size)
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

# Z4 = CyclicGroup(4)
# e = GroupElement(0, Z4)
# a = GroupElement(1, Z4)

# A = Obj(e=>2, a=>2, a*a=>3, a*a*a=>4)
# B = Obj(e=>2, a=>4, a*a=>3, a*a*a=>2)

# A_dual = dual_obj(A)
# B_dual = dual_obj(B)

# T1 = random_mor(Float64, (A,B,A,B))
# T2 = random_mor(Float64, (B_dual, A_dual, A,B))

# T = VecG_tensordot(T1, T2, 2)
include("groups.jl")
include("block_matrix_calculation.jl")
include("VecGtensor.jl")
include("display.jl")

function extract_blocks_to_matrix(mor::Mor{G, T}, n_leg_split::Int, g_bridge::GroupElement{G}) where {T, G <: Group} # legs split as 1, ..., n-1, n, |, n+1, n+2, ...
    n_leg = length(mor.objects)
    group = get_group(mor)

    # Check if the elements of `legs` are within the range 1 to n
    if n_leg_split<1 || n_leg_split>n_leg
        error("Elements of legs must be in the range 1 to $n_leg")
    end

    out_sectors = group_tree(g_bridge, n_leg_split)
    in_sectors = group_tree(inverse(g_bridge), n_leg - n_leg_split)

    A_blocks = Matrix{AbstractMatrix}(undef, length(out_sectors), length(in_sectors))

    for (i,out_sector) in enumerate(out_sectors), (j,in_sector) in enumerate(in_sectors)
        tot_sect = Sector(out_sector..., in_sector...)
        # @show tot_sect
        dims = get_sector_size(mor, tot_sect)
        row_dims = prod(dims[1:n_leg_split])
        col_dims = prod(dims[n_leg_split+1:n_leg])
        # @show dims, row_dims, col_dims
        if haskey(mor.data, tot_sect)
            tensor = mor.data[tot_sect]
            matrix = reshape(tensor, row_dims, col_dims)
            A_blocks[i, j] = matrix
        else
            A_blocks[i, j] = zeros(row_dims, col_dims)
        end
    end
    return A_blocks
end

function VecG_factorize(mor::Mor{G, T}, n_leg_split::Int, Dcut::Int, method::AbstractString) where {T, G<:Group}
    group = get_group(mor)
    n_leg = length(mor.objects)
    Fobj = (mor[1:n_leg_split]...,zero_obj(group))
    Kobj = (mor[n_leg_split+1:end]...,zero_obj(group))
    F = Mor(T, Fobj)
    K = Mor(T, Kobj)
    for g_bridge in elements(group)
        block_matrix = extract_blocks_to_matrix(mor, n_leg_split, g_bridge)
        F_mat, K_mat = factorize_block_matrix(block_matrix, method, Dcut)
        if F_mat == undef || K_mat == undef
            F[end][inverse(g_bridge)] = 0
            K[end][g_bridge] = 0
            out_sectors = group_tree(g_bridge, n_leg_split)
            in_sectors = group_tree(inverse(g_bridge), n_leg - n_leg_split)
            for out_sector in out_sectors
                F[out_sector..., inverse(g_bridge)] = zeros(get_sector_size(F, (out_sector..., inverse(g_bridge)))...)
            end
            for in_sector in in_sectors
                K[in_sector..., g_bridge] = zeros(get_sector_size(K, (in_sector..., g_bridge))...)
            end
        else
            multiplicity = size(F_mat[1],2)
            F[end][inverse(g_bridge)] = multiplicity
            K[end][g_bridge] = multiplicity
            out_sectors = group_tree(g_bridge, n_leg_split)
            in_sectors = group_tree(inverse(g_bridge), n_leg - n_leg_split)
            for (i,out_sector) in enumerate(out_sectors)
                F[out_sector..., inverse(g_bridge)] = reshape(F_mat[i], get_sector_size(F, (out_sector..., inverse(g_bridge))))
            end
            for (j,in_sector) in enumerate(in_sectors)
                K[in_sector..., g_bridge] = reshape(K_mat[j], get_sector_size(K, (in_sector..., g_bridge)))
            end
        end
    end
    return F, K
end




function is_accend(tup::Tuple{Vararg{Int}}, modn::Int)
    test = true
    for i = 1 : length(tup)-1
        if !(tup[i] in 1:modn) || !(tup[i+1] in 1:modn) || tup[i+1] != (mod(tup[i], modn) + 1)
            test = false
        end
    end
    return test
end

function is_cyclic(tup::Tuple{Vararg{Int}}, modn::Int)
    test = is_accend(tup, modn)
    if tup[1]!=(mod(tup[end], modn) + 1)
        test = false
    end
    return test
end

function to_perm(tup::Tuple{Vararg{Int}}, modn::Int)
    if !is_accend(tup, modn)
        throw(Argumenterror("The $tup is not in an accending order"))
    end
    newtup = tup[1]:tup[1]+modn-1
    mod_newtup = Tuple(map(x->mod(x-1,modn)+1, newtup))
    return mod_newtup
end

function VecG_permutesectors(sect::Sector, perm::Tuple{Vararg{Int}})
    modn = length(sect.sect)
    if is_cyclic(perm, modn) == false
        throw(ArgumentError("The permutation $perm is not cyclic"))
    end
    perm_sect_tup = ()
    for i in perm
        perm_sect_tup = (perm_sect_tup..., sect[i])
    end
    return Sector(perm_sect_tup...)
end

function VecG_permutedims(mor::Mor{G, T}, perm::Tuple{Vararg{Int}}) where {T, G<:Group}
    modn = length(mor.objects)
    if is_cyclic(perm, modn) == false
        throw(ArgumentError("The permutation $perm is not cyclic"))
    end
    perm_obj = ()
    for perm_i in perm
        perm_obj = (perm_obj..., mor[perm_i])
    end
    perm_mor = Mor(T, perm_obj)
    for sect in keys(mor.data)
        perm_sect = VecG_permutesectors(sect, perm)
        perm_tensor = permutedims(mor[sect], perm)
        perm_mor[perm_sect] = perm_tensor
    end
    return perm_mor
end

function VecG_factorize(mor::Mor, n_leg_split::Tuple{Vararg{Int}}, Dcut::Int, method::AbstractString)
    modn = length(mor.objects)
    if is_accend(n_leg_split, modn) == false
        throw(ArgumentError("The factorize leg $n_leg_split is not accending"))
    end

    perm = to_perm(n_leg_split, modn)
    perm_mor = VecG_permutedims(mor, perm)

    F, K = VecG_factorize(perm_mor, length(n_leg_split), Dcut, method)

    return F, K
end


function VecG_tensordot(A::Mor{G,T}, B::Mor{G,T}, leg_cont::Int) where {T, G<:Group}
    A_legs = length(A.objects)
    B_legs = length(B.objects)

    for i in 1:n
        if A[end-i+1]!=B[i]
            throw(ArgumentError("Contract index is illegal"))
        end
    end

    newobjs = (A[1:A_legs -leg_cont]..., B[leg_cont+1:B_legs])

    Cont = Mor(T, newobjs)
end

D4 = DihedralGroup(4)
s = GroupElement((1,0),D4)
r = GroupElement((0,1),D4)
e = identity(D4)
A = Obj(e=>2, s=>2, r=>3)
B = Obj(e=>2, s*r=>3, s=>4)
# get_group(A)
# zero_obj(D4)

T = random_mor(Float64, (A, A, B, B))
# F, K = VecG_factorize(T, 2, 10, "qr")

# sect = Sector(e, s, r, s*r)
# VecGpermutesectors(sect, (3,4,1,2))

# VecGpermutedims(T, (2,3,4,1))

F1, K1 = VecG_factorize(T, (1,2), 20, "svd")
F2, K2 = VecG_factorize(T, (2,3), 20, "svd")
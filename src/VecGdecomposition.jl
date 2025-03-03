"""
This module implements decompositions for VecG-tensors.
"""

"""
Given a bridge group element g_bridge, find all sectors that fuses to g_bridge when splited to 1, ..., n_leg_split-1, n_leg_split, | n_leg_split+1, n_leg_split+2, ...

Graphically, we consider the following:

```julia
|  |  |  |  |
g1 g2...    g_n_leg_split
|  |  |  |  |
^  ^  ^  ^  ^
|  |  |  |  |
 \\  \\ | /  /
      ^
      |
      g_bridge
      |
 /  / | \\  \\
|  |  |  |  |
gn gn-1...g_n_leg_split+1  
|  |  |  |  |
v  v  v  v  v
|  |  |  |  |
```

# Input:
- a morphism
- an integer that describes the number of splitted legs

# Output:
- a dictionary, whose keys are bridge sectors appeared in the morphism, and values are a 2-argument tuple, each of which is also a tuple, representing the out sectors and in sectors

# Example

```
    D6 = DihedralGroup(3)
    e = identity_element(D6)
    s = GroupElement((1,0), D6)
    r = GroupElement((0,1), D6)
    A = Obj(e=>2, s=>3, r=>2, s*r=>15) 
    B = Obj(e=>2, s*r*r=>4, r=>3, s*r=>2)
    C = Obj(e=>2, s=>3, r*r=>2, s*r=>2)
    D = Obj(e=>2, s=>4, r=>3, s*r*r=>15)
    TT = random_mor(Float64, (A, B, C, D))
    to_sector_outin(TT, 2)
    Dict{Any, Any} with 6 entries:
    s   => (((s, e), (r, sr)), ((r², sr²), (s, e), (e, s)))
    e   => (((e, e), (sr, sr)), ((e, e), (r², r), (s, s)))
    sr  => (((r, sr²), (sr, e), (s, r), (e, sr)), ((sr, e), (s, r), (r², s)))
    r   => (((r, e), (sr, sr²), (e, r), (s, sr)), ((s, sr²), (r², e), (sr, s)))
    r²  => (((s, sr²), (r, r)), ((e, r), (sr, sr²)))
    sr² => (((sr, r), (e, sr²)), ((e, sr²), (sr, r)))
```
"""
function to_sector_outin(mor::Mor{G, T}, n_leg_split::Int) where {T, G <: Group}
    sector_outin = Dict()
    for sector in keys(mor.data)
        out_sect = sector.sect[1:n_leg_split]
        in_sect = sector.sect[n_leg_split+1:end]
        g_bridge = multiply(out_sect)
        if haskey(sector_outin, g_bridge)
            out_sector_tuple, in_sector_tuple = sector_outin[g_bridge]
            if !(out_sect in out_sector_tuple)
                out_sector_tuple = (out_sector_tuple..., out_sect)
            end
            if !(in_sect in in_sector_tuple)
                in_sector_tuple = (in_sector_tuple..., in_sect)
            end
            sector_outin[g_bridge] = (out_sector_tuple, in_sector_tuple)
        else
            sector_outin[g_bridge] = ((out_sect,), (in_sect,))
        end
    end

    return sector_outin
end


"""
Extracts VecG morphism to block matrices when given a fixed group element g_bridge. This is an intermediate function before performing svd.

# Input
- A morphism
- A integer, before and include which legs will be considered as row indices, and after which legs will be considered as column indices
- A dictionary that representing the splitting of the sector, produced from to_sector_outin

# Output:
- A dictionary, whose keys are bridge sectors appeared in the morphism, and values are matrix of matrices

# Example

```julia
julia> D6 = DihedralGroup(3)
        e = identity_element(D6)
        s = GroupElement((1,0), D6)
        r = GroupElement((0,1), D6)
        A = Obj(e=>2, s=>3, r=>2, s*r=>1) 
        B = Obj(e=>2, s*r*r=>2, r=>3, s*r=>2)
        C = Obj(e=>2, s=>3, r*r=>2, s*r=>2)
        D = Obj(e=>2, s=>4, r=>1, s*r*r=>1)
        TT = random_mor(Float64, (A, B, C, D))
        sec_split = to_sector_outin(TT, 2)
        Mat = extract_blocks_to_matrix(TT, 2, sec_split)
        @show Mat[s][1,1] == reshape(TT[s,e,r*r,s*r*r], 6,2) # true
        @show Mat[s][1,2] == reshape(TT[s,e,s,e], 6,6) # true
        @show Mat[s][2,2] == reshape(TT[r,s*r,s,e], 4,6) # true

```
"""
function extract_blocks_to_matrix(mor::Mor{G, T}, n_leg_split::Int, sector_outin) where {T, G <: Group} # legs split as 1, ..., n-1, n, |, n+1, n+2, ...
    n_leg = length(mor.objects)
    block_matrix_dic = Dict()

    for g in keys(sector_outin)
        out_sector_tuple, in_sector_tuple = sector_outin[g]

        A_blocks = Matrix{AbstractMatrix}(undef, length(out_sector_tuple), length(in_sector_tuple))

        for (i,out_sector) in enumerate(out_sector_tuple), (j,in_sector) in enumerate(in_sector_tuple)
            tot_sect = Sector(out_sector..., in_sector...)
            dims = get_sector_size(mor, tot_sect)
            row_dims = prod(dims[1:n_leg_split])
            col_dims = prod(dims[n_leg_split+1:n_leg])
            if haskey(mor.data, tot_sect)
                tensor = mor.data[tot_sect]
                matrix = reshape(tensor, row_dims, col_dims)
                A_blocks[i, j] = matrix
            else
                A_blocks[i, j] = zeros(row_dims, col_dims)
            end
        end
        block_matrix_dic[g] = A_blocks
    end
    
    return block_matrix_dic
end

function svd_calculation(block_matrix_dic)
    svd_dic = Dict()
    for g in keys(block_matrix_dic)
        svd_dic[g] = block_matrix_svd(block_matrix_dic[g])
    end
    return svd_dic
end

function merge_matrix_to_block(mor::Mor{G, T}, n_leg_split::Int, svd_dic, sector_outin) where {T, G<:Group}
    svd_leg_V = Dict{GroupElement{G}, Int}()
    for g in keys(svd_dic)
        svd_leg_V[g] = length(svd_dic[g][2])
    end
    obj_bridge_V = Obj(svd_leg_V)
    obj_bridge_U = dual_obj(obj_bridge_V)

    U = Mor(T, (mor[1:n_leg_split]...,obj_bridge_U))
    V = Mor(T, (mor[n_leg_split+1:end]...,obj_bridge_V))
    S = Mor(Float64, (obj_bridge_V, obj_bridge_U))

    for g_bridge in keys(svd_dic)
        U_mat, S_vec, V_mat = svd_dic[g_bridge]
        out_sector_tuple, in_sector_tuple = sector_outin[g_bridge]
        for (i,out_sector) in enumerate(out_sector_tuple)
            U[out_sector..., inverse(g_bridge)] = reshape(U_mat[i], get_sector_size(U, (out_sector..., inverse(g_bridge))))
        end
        for (j,in_sector) in enumerate(in_sector_tuple)
            V[in_sector..., g_bridge] = reshape(conj.(V_mat[j]), get_sector_size(V, (in_sector..., g_bridge))) # Here conjugation is taken to match the tensor network convension
        end
        S[g_bridge,inverse(g_bridge)] = diagm(S_vec)
    end

    return U, S, V
end

function block_cutoff!(svd_dic, Dcut::Int)
    S_tot = []
    
    for g in keys(svd_dic)
        append!(S_tot, svd_dic[g][2])
    end

    sorted_singular_values = sort(abs.(S_tot), rev=true)  # 从大到小排序
    cutoff = min(Dcut, length(S_tot))
    cutoff_threshold = sorted_singular_values[cutoff]

    for g in keys(svd_dic)
        U, S, V = svd_dic[g]
        Dcut_sect = sum(abs.(S) .>= cutoff_threshold)
        newS = S[1:Dcut_sect]
        newU = []
        newV = []
        for u in U
            push!(newU, u[:,1:Dcut_sect])
        end
        for v in V
            push!(newV, v[:,1:Dcut_sect])
        end
        svd_dic[g] = (newU, newS, newV)
    end
    
    return svd_dic
end

function VecG_svd(mor::Mor, n_leg_split::Tuple{Vararg{Int}}, Dcut::Int)
    modn = length(mor.objects)
    split_number = length(n_leg_split)
    if is_accend(n_leg_split, modn) == false
        throw(ArgumentError("The factorize leg $n_leg_split is not accending"))
    end

    perm = to_perm(n_leg_split, modn)
    perm_mor = VecG_permutedims(mor, perm)

    sec_split = to_sector_outin(perm_mor, split_number)
    Mat = extract_blocks_to_matrix(perm_mor, split_number, sec_split)
    svd_dic = svd_calculation(Mat)
    block_cutoff!(svd_dic, Dcut)
    U, S, V = merge_matrix_to_block(perm_mor, split_number, svd_dic, sec_split)

    return U, S, V
end

function VecG_qr(mor::Mor, n_leg_split::Tuple{Vararg{Int}})
    modn = length(mor.objects)
    if is_accend(n_leg_split, modn) == false
        throw(ArgumentError("The factorize leg $n_leg_split is not accending"))
    end

    perm = to_perm(n_leg_split, modn)
    perm_mor = VecG_permutedims(mor, perm)

    Q, R = VecG_qr(perm_mor, length(n_leg_split))

    return Q, R
end

function VecG_qr(mor::Mor{G, T}, n_leg_split::Int) where {T, G<:Group}
    group = get_group(mor)
    n_leg = length(mor.objects)
    Q = Mor(T, (mor[1:n_leg_split]...,zero_obj(group)))
    R = Mor(T, (zero_obj(group), mor[n_leg_split+1:end]...))

    for g_bridge in elements(group)
        block_matrix = extract_blocks_to_matrix(mor, n_leg_split, g_bridge)
        if 0 in size(block_matrix)
            throw(ArgumentError("The QR decomposition cannot be done if the matrix has zero column"))
        else
            Q_mat, R_mat = block_matrix_qr(block_matrix)
            multiplicity = size(Q_mat[1], 2)
            Q[end][inverse(g_bridge)] = multiplicity
            R[1][g_bridge] = multiplicity
            out_sectors = group_tree(g_bridge, n_leg_split)
            in_sectors = group_tree(inverse(g_bridge), n_leg - n_leg_split)
            for (i,out_sector) in enumerate(out_sectors)
                Q[out_sector..., inverse(g_bridge)] = reshape(Q_mat[i], get_sector_size(Q, (out_sector..., inverse(g_bridge))))
            end
            for (j,in_sector) in enumerate(in_sectors)
                R[g_bridge, in_sector...] = reshape(R_mat[j], get_sector_size(R, (g_bridge, in_sector...)))
            end
        end
    end
    return Q, R
end
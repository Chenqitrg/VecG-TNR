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
- a KeySet
- an integer that describes the number of splitted legs
- a group element that describes the total sector of splitting

# Output:
- a tuple of group element tuples, describing the out sector
- a tuple of group element tuples, describing the in sector

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
    @show to_sector_outin(keys(TT.data), 2, s)
    (((s, e), (r, sr)), ((e, s), (s, e), (r², sr²)))
```
"""
function to_sector_outin(ks::Base.KeySet, n_leg_split::Int, g_bridge::GroupElement{G}) where {G <: Group}
    out_sector = ()
    in_sector = ()

    for sector in ks
        out_sect = sector.sect[1:n_leg_split]
        in_sect = sector.sect[n_leg_split+1:end]
        if multiply(out_sect) == g_bridge
            if !(out_sect in out_sector)
                out_sector = (out_sector..., out_sect)
            end
            if !(in_sect in in_sector)
                in_sector = (in_sector..., in_sect)
            end
        end
    end

    return out_sector, in_sector
end



"""
Extracts VecG morphism to block matrices when given a fixed group element g_bridge. This is an intermediate function before performing svd.

# Input
- A morphism
- A integer, before and include which legs will be considered as row indices, and after which legs will be considered as column indices
- A group element, representing the total sector to be split. Only the sectors such that rows are multiplied to that group element is taken out.

# Output:
- A matrix of matrices


"""
function extract_blocks_to_matrix(mor::Mor{G, T}, out_sector, in_sector) where {T, G <: Group} # legs split as 1, ..., n-1, n, |, n+1, n+2, ...
    n_leg = length(mor.objects)
    group = get_group(mor)

    # Check if the elements of `legs` are within the range 1 to n
    if n_leg_split<1 || n_leg_split>n_leg
        error("Elements of legs must be in the range 1 to $n_leg")
    end

   
    # out_sectors = group_tree(g_bridge, n_leg_split)
    # in_sectors = group_tree(inverse(g_bridge), n_leg - n_leg_split)

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

function VecG_svd(mor::Mor{G, T}, n_leg_split::Int) where {T, G<:Group}
    group = get_group(mor)
    n_leg = length(mor.objects)
    U = Mor(T, (mor[1:n_leg_split]...,zero_obj(group)))
    V = Mor(T, (mor[n_leg_split+1:end]...,zero_obj(group)))
    S = Mor(T, (zero_obj(group),zero_obj(group)))
    for g_bridge in elements(group)
        block_matrix = extract_blocks_to_matrix(mor, n_leg_split, g_bridge)
        U_mat, S_vec, V_mat = block_matrix_svd(block_matrix)
        multiplicity = length(S_vec)
        U[end][inverse(g_bridge)] = multiplicity
        V[end][g_bridge] = multiplicity
        S[1][g_bridge] = multiplicity
        S[2][inverse(g_bridge)] = multiplicity
        out_sectors = group_tree(g_bridge, n_leg_split)
        in_sectors = group_tree(inverse(g_bridge), n_leg - n_leg_split)
        for (i,out_sector) in enumerate(out_sectors)
            U[out_sector..., inverse(g_bridge)] = reshape(U_mat[i], get_sector_size(U, (out_sector..., inverse(g_bridge))))
        end
        for (j,in_sector) in enumerate(in_sectors)
            V[in_sector..., g_bridge] = reshape(conj.(V_mat[j]), get_sector_size(V, (in_sector..., g_bridge)))
        end
        S[g_bridge,inverse(g_bridge)] = T.(diagm(S_vec))
    end

    return U, S, V

end

function VecG_cutoff(U::Mor{G,T}, S::Mor{G,T}, V::Mor{G,T}, Dcut::Int) where {T,G<:Group}
    group = get_group(S)
    S_tot = T[]
    out_legs = length(U.objects) - 1
    in_legs = length(V.objects) - 1
    
    for g in elements(group)
        append!(S_tot, diag(S[g, inverse(g)]))
    end

    S_tot_vcat = vcat(S_tot)
    
    sorted_singular_values = sort(abs.(S_tot), rev=true)  # 从大到小排序
    cutoff = min(Dcut, length(S_tot_vcat))
    cutoff_threshold = sorted_singular_values[cutoff]

    for g in elements(group)
        for out_sect in group_tree(g, out_legs), in_sect in group_tree(inverse(g), in_legs)
            Dcut_sect = sum(abs.(diag(S[g, inverse(g)])) .>= cutoff_threshold)
            U_indices = ntuple(_ -> :, out_legs)
            V_indices = ntuple(_->:, in_legs)
            U[end][inverse(g)] = Dcut_sect
            V[end][g] = Dcut_sect
            S[1][g] = Dcut_sect
            S[2][inverse(g)] = Dcut_sect
            U[out_sect..., inverse(g)] = U[out_sect..., inverse(g)][U_indices..., 1:Dcut_sect]
            V[in_sect..., g] = V[in_sect..., g][V_indices..., 1:Dcut_sect]
            S[g, inverse(g)] = S[g, inverse(g)][1:Dcut_sect, 1:Dcut_sect]
        end
    end
    
    return U, S, V
end

function VecG_cutoff(U::Mor{G,T}, S::Mor{G,T}, V::Mor{G,T}, epsilon::Float64) where {T,G<:Group}
    group = get_group(S)

    out_legs = length(U.objects) - 1
    in_legs = length(V.objects) - 1

    for g in elements(group)
        for out_sect in group_tree(g, out_legs), in_sect in group_tree(inverse(g), in_legs)
            Dcut_sect = sum(abs.(diag(S[g, inverse(g)])) .>= epsilon)
            U_indices = ntuple(_ -> :, out_legs)
            V_indices = ntuple(_->:, in_legs)
            U[end][inverse(g)] = Dcut_sect
            V[end][g] = Dcut_sect
            S[1][g] = Dcut_sect
            S[2][inverse(g)] = Dcut_sect
            U[out_sect..., inverse(g)] = U[out_sect..., inverse(g)][U_indices..., 1:Dcut_sect]
            V[in_sect..., g] = V[in_sect..., g][V_indices..., 1:Dcut_sect]
            S[g, inverse(g)] = S[g, inverse(g)][1:Dcut_sect, 1:Dcut_sect]
        end
    end
    
    return U, S, V
end

function VecG_cutoff(U::Mor{G,T}, S::Mor{G,T}, V::Mor{G,T}, Dcut::Int, epsilon::Float64) where {T,G<:Group}
    group = get_group(S)
    S_tot = T[]
    out_legs = length(U.objects) - 1
    in_legs = length(V.objects) - 1
    
    for g in elements(group)
        append!(S_tot, diag(S[g, inverse(g)]))
    end

    S_tot_vcat = vcat(S_tot)
    
    sorted_singular_values = sort(abs.(S_tot), rev=true)  # 从大到小排序
    cutoff = min(Dcut, length(S_tot_vcat))
    cutoff_threshold = max(sorted_singular_values[cutoff],epsilon)

    for g in elements(group)
        for out_sect in group_tree(g, out_legs), in_sect in group_tree(inverse(g), in_legs)
            Dcut_sect = sum(abs.(diag(S[g, inverse(g)])) .>= cutoff_threshold)
            U_indices = ntuple(_ -> :, out_legs)
            V_indices = ntuple(_->:, in_legs)
            U[end][inverse(g)] = Dcut_sect
            V[end][g] = Dcut_sect
            S[1][g] = Dcut_sect
            S[2][inverse(g)] = Dcut_sect
            U[out_sect..., inverse(g)] = U[out_sect..., inverse(g)][U_indices..., 1:Dcut_sect]
            V[in_sect..., g] = V[in_sect..., g][V_indices..., 1:Dcut_sect]
            S[g, inverse(g)] = S[g, inverse(g)][1:Dcut_sect, 1:Dcut_sect]
        end
    end
    
    return U, S, V
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

function VecG_svd(mor::Mor, n_leg_split::Tuple{Vararg{Int}}, Dcut::Int)
    modn = length(mor.objects)
    if is_accend(n_leg_split, modn) == false
        throw(ArgumentError("The factorize leg $n_leg_split is not accending"))
    end

    perm = to_perm(n_leg_split, modn)
    perm_mor = VecG_permutedims(mor, perm)

    U, S, V = VecG_svd(perm_mor, length(n_leg_split))
    U, S, V = VecG_cutoff(U, S, V, Dcut)

    return U, S, V
end

function VecG_svd(mor::Mor, n_leg_split::Tuple{Vararg{Int}}, epsilon::Float64)
    modn = length(mor.objects)
    if is_accend(n_leg_split, modn) == false
        throw(ArgumentError("The factorize leg $n_leg_split is not accending"))
    end

    perm = to_perm(n_leg_split, modn)
    perm_mor = VecG_permutedims(mor, perm)

    U, S, V = VecG_svd(perm_mor, length(n_leg_split))
    U, S, V = VecG_cutoff(U, S, V, epsilon)

    return U, S, V
end

function VecG_svd(mor::Mor, n_leg_split::Tuple{Vararg{Int}}, Dcut::Int, epsilon::Float64)
    modn = length(mor.objects)
    if is_accend(n_leg_split, modn) == false
        throw(ArgumentError("The factorize leg $n_leg_split is not accending"))
    end

    perm = to_perm(n_leg_split, modn)
    perm_mor = VecG_permutedims(mor, perm)

    U, S, V = VecG_svd(perm_mor, length(n_leg_split))
    U, S, V = VecG_cutoff(U, S, V, Dcut, epsilon)

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

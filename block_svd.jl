include("groups.jl")
include("block_matrix_calculation.jl")
include("VecGtensor.jl")


function extract_blocks_to_matrix(bt::BlockTensor{T, G}, n_leg_split::Int, g_bridge::GroupElement{G}) where {T, G <: Group} # legs split as 1, ..., n-1, n, |, n+1, n+2, ...
    n_leg = length(bt.objects)
    group = bt.group

    # Check if the elements of `legs` are within the range 1 to n
    if n_leg_split<1 || n_leg_split>n_leg
        error("Elements of legs must be in the range 1 to $n_leg")
    end

    out_sectors = group_tree(group, g_bridge, n_leg_split)
    in_sectors = group_tree(group, inverse(g_bridge), n_leg - n_leg_split)

    A_blocks = Matrix{Any}(undef, length(out_sectors), length(in_sectors))

    for (i,out_sector) in enumerate(out_sectors), (j,in_sector) in enumerate(in_sectors)
        tot_sect = (out_sector..., in_sector...)
        @show tot_sect
        dims = sector_size(bt, tot_sect)
        row_dims = prod(dims[1:n_leg_split])
        col_dims = prod(dims[n_leg_split+1:n_leg])
        @show dims, row_dims, col_dims
        if haskey(bt.data, tot_sect)
            tensor = bt.data[tot_sect]
            matrix = reshape(tensor, row_dims, col_dims)
            A_blocks[i, j] = matrix
        else
            A_blocks[i, j] = zeros(row_dims, col_dims)
        end
    end
    return A_blocks
end

function block_factorize(bt::BlockTensor{T, G}, n_leg_split::Int, Dcut::Int, method::AbstractString) where {T, G <: Group}
    group = bt.group
    
end


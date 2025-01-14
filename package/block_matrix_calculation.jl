function concatenate_matrices_with_metadata(matrices::AbstractMatrix{<:AbstractMatrix})
    # Number of submatrices in rows and columns
    rows::Int = size(matrices, 1)  # Number of submatrices in each row
    cols::Int = size(matrices, 2)  # Number of submatrices in each column
    
    # Dictionary to store metadata: (i, j) => (start_row, start_col, num_rows, num_cols)
    metadata::Dict{Tuple{Int, Int}, Tuple{Int, Int, Int, Int}} = Dict()
    
    # Compute row and column sizes of the final large matrix
    row_sizes::Vector{Int} = [maximum(size(matrices[i, j], 1) for j in 1:cols) for i in 1:rows]
    col_sizes::Vector{Int} = [maximum(size(matrices[i, j], 2) for i in 1:rows) for j in 1:cols]
    
    total_rows::Int = sum(row_sizes)  # Total number of rows in the large matrix
    total_cols::Int = sum(col_sizes)  # Total number of columns in the large matrix
    
    # Create the final large matrix
    big_matrix = Matrix{eltype(matrices[1, 1])}(undef, total_rows, total_cols)
    
    # Fill the large matrix and record metadata
    current_row::Int = 1
    for i in 1:rows
        current_col::Int = 1
        for j in 1:cols
            sub_matrix = matrices[i, j]
            num_rows, num_cols = size(sub_matrix)
            
            # Insert the submatrix into the large matrix
            big_matrix[current_row:current_row+num_rows-1, current_col:current_col+num_cols-1] .= sub_matrix
            
            # Record metadata for the submatrix
            metadata[(i, j)] = (current_row, current_col, num_rows, num_cols)
            
            # Update the starting column for the next submatrix
            current_col += col_sizes[j]
        end
        # Update the starting row for the next row of submatrices
        current_row += row_sizes[i]
    end
    
    return big_matrix, metadata
end


function block_matrix_svd(matrices::AbstractMatrix{<:AbstractMatrix})
    # Extract row sizes and column sizes for the block structure
    r_sizes = [size(matrices[i, 1], 1) for i in 1:size(matrices, 1)]  # Row sizes of each block row
    c_sizes = [size(matrices[1, j], 2) for j in 1:size(matrices, 2)]  # Column sizes of each block column

    # Concatenate input block matrices into a single large matrix
    matrices_big, _ = concatenate_matrices_with_metadata(matrices)

    # Perform matrix factorization

    U, S, V = svd(matrices_big)
    # Split F into block rows according to r_sizes
    U_blocks = Vector{AbstractMatrix}(undef, length(r_sizes))
    row_start = 1
    for (i, r) in enumerate(r_sizes)
        U_blocks[i] = U[row_start:row_start + r - 1, :]
        row_start += r
    end
    # Split K into block columns according to c_sizes
    V_blocks = Vector{AbstractMatrix}(undef, length(c_sizes))
    col_start = 1
    for (j, c) in enumerate(c_sizes)
        V_blocks[j] = V[col_start:col_start + c - 1, :]
        col_start += c
    end
    return U_blocks, S, V_blocks
end

function block_matrix_qr(matrices::AbstractMatrix{<:AbstractMatrix})
    # Extract row sizes and column sizes for the block structure
    r_sizes = [size(matrices[i, 1], 1) for i in 1:size(matrices, 1)]  # Row sizes of each block row
    c_sizes = [size(matrices[1, j], 2) for j in 1:size(matrices, 2)]  # Column sizes of each block column

    # Concatenate input block matrices into a single large matrix
    matrices_big, _ = concatenate_matrices_with_metadata(matrices)
    Q, R = qr(matrices_big)
    Rrow = size(R,1)
    Q = Q[:,1:Rrow]
    # Split F into block rows according to r_sizes
    Q_blocks = Vector{AbstractMatrix}(undef, length(r_sizes))
    row_start = 1
    for (i, r) in enumerate(r_sizes)
        Q_blocks[i] = Q[row_start:row_start + r - 1, :]
        row_start += r
    end

    # Split K into block columns according to c_sizes
    R_blocks = Vector{AbstractMatrix}(undef, length(c_sizes))
    col_start = 1
    for (j, c) in enumerate(c_sizes)
        R_blocks[j] = R[:,col_start:col_start + c - 1]
        col_start += c
    end

    # Return the block matrices Q and R
    return Q_blocks, R_blocks
end

function test_block_matrix_svd()
    # Define submatrices
    A = [1 2; 3 4]  # 2x3 matrix
    B = [5 6 6; 7 8 8]  # 2x4 matrix
    C = [9 10; 11 12; 15 16; 17 18]  # 3x3 matrix
    D = [13 14 19; 15 16 20; 21 22 23; 25 26 27]  # 3x4 matrix

    # Create a 2x2 "matrix of matrices" explicitly
    matrices = Array{Matrix{Float64}}(undef, 2, 2)  # Create a 2x2 array to hold submatrices
    matrices[1, 1] = A
    matrices[1, 2] = B
    matrices[2, 1] = C
    matrices[2, 2] = D

    # # Call the function
    # big_matrix, metadata = concatenate_matrices_with_metadata(matrices)
    return block_matrix_svd(matrices)
end

function partial_trace(arr::Array{T}, leg_cont::Int) where {T}
    sizeT = size(arr)
    tot_legs = length(sizeT)
    num_remain = tot_legs - 2 * leg_cont
    new_arr = zeros(T, sizeT[1:num_remain])
    new_iter = Iterators.product((1:n for n in sizeT[1:num_remain])...)
    sum_iter = Iterators.product((1:n for n in sizeT[num_remain+1:num_remain+leg_cont])...)
    for newind in new_iter
        for reind in sum_iter
            new_arr[newind...] = new_arr[newind...] + arr[newind..., reind..., reverse(reind)...]
        end
    end
    return new_arr
end

function test_block_matrix_qr()
    # Define submatrices
    A = [1 2; 3 4]  # 2x3 matrix
    B = [5 6 6; 7 8 8]  # 2x4 matrix
    C = [9 10; 11 12; 15 16; 17 18]  # 3x3 matrix
    D = [13 14 19; 15 16 20; 21 22 23; 25 26 27]  # 3x4 matrix

    # Create a 2x2 "matrix of matrices" explicitly
    matrices = Array{Matrix{Float64}}(undef, 2, 2)  # Create a 2x2 array to hold submatrices
    matrices[1, 1] = A
    matrices[1, 2] = B
    matrices[2, 1] = C
    matrices[2, 2] = D

    # # Call the function
    # big_matrix, metadata = concatenate_matrices_with_metadata(matrices)
    return block_matrix_qr(matrices)
end

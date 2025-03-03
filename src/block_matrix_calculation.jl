
"""
To perform the calculation with VecG-tensor, we need to reduce the formal calculation to the detailed matrix calculation. 
This module provides functions for block matrix calculations, including block matrix SVD and QR decomposition.
    It also provides functions for partial trace and outer product of tensors.
        Our philosophy is to separate the formal part of the code from the numerical part.
        The formal part is in the VecG_TNR module, and the numerical part is in the main module.
        The formal part is for the formal definition of the tensor network renormalization, which is mostly the structure of the data.
        And the numerical part is for the numerical implementation, which is mostly the calculation of the data.
"""

"""
Concatenate a 2D array of matrices into a single large matrix, and record metadata for the submatrices.
    The metadata is stored in a dictionary with keys (i, j), which represents the location of the block.
    And values (start_row, start_col, num_rows, num_cols) represent the starting row and column of the block, and the number of rows and columns of the block.
    The metadata can be used to extract the submatrices from the large matrix.

# Arguments
- `matrices::AbstractMatrix{<:AbstractMatrix}`: A 2D array of matrices to be concatenated.

# Returns
- `big_matrix::Matrix`: The large matrix obtained by concatenating the input matrices.
- `metadata::Dict{Tuple{Int, Int}, Tuple{Int, Int, Int, Int}}`: The metadata for the submatrices.

# Example
```julia
M = Matrix{Matrix}(undef, 2, 2)
A = [1 2; 3 4]  # 2x2 matrix
B = [5 6 6; 7 8 8]  # 2x3 matrix
C = [9 10; 11 12; 15 16; 17 18]  # 3x2 matrix
D = [13 14 19; 15 16 20; 21 22 23]  # 3x3 matrix
M[1,1] = A
M[1,2] = B
M[2,1] = C
M[2,2] = D
big_matrix, metadata = concatenate_matrices_with_metadata(M)
@show big_matrix
6×5 Matrix{Int64}:
  1   2   5   6   6
  3   4   7   8   8
  9  10  13  14  19
 11  12  15  16  20
 15  16  21  22  23
 17  18   0   0   0
 @show metadata
 Dict{Tuple{Int64, Int64}, NTuple{4, Int64}} with 4 entries:
  (1, 2) => (1, 3, 2, 3)
  (1, 1) => (1, 1, 2, 2)
  (2, 2) => (3, 3, 3, 3)
  (2, 1) => (3, 1, 4, 2)
```
"""
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

"""
Perform SVD on a block matrix, where each block is a submatrix.
    The input is a 2D array of matrices, and the output is a tuple of block matrices U, S, and V.
    The block matrices U and V are the left and right singular vectors, respectively.
    The block matrix S is the singular values, stored as a vector.

# Arguments
- `matrices::AbstractMatrix{<:AbstractMatrix}`: A 2D array of matrices to be factorized.

# Returns
- `U_blocks::Vector{AbstractMatrix}`: A vector of block matrices representing the left singular vectors.
- `S::Vector`: A vector of singular values.
- `V_blocks::Vector{AbstractMatrix}`: A vector of block matrices representing the right singular vectors.

# Example
```julia
A = [1 2; 3 4]  # 2x3 matrix
B = [5 6 6; 7 8 8]  # 2x3 matrix
C = [9 10; 11 12; 15 16; 17 18]  # 4x2 matrix
D = [13 14 19; 15 16 20; 21 22 23; 25 26 27]  # 4x3 matrix
matrices = Array{Matrix{Float64}}(undef, 2, 2)
matrices[1, 1] = A
matrices[1, 2] = B
matrices[2, 1] = C
matrices[2, 2] = D
U_blocks, S, V_blocks = block_matrix_svd(matrices)
```

The output U_blocks, S, and V_blocks can be compared with the output of the svd function. Here we can use the concatenate_matrices_with_metadata function to concatenate the block matrices into a single large matrix and compare the results.
However, before doing so, we need to reshape the block matrices into a 2x1 array to match the input format of the concatenate_matrices_with_metadata function.
The comparison should return true for both U and V.

```julia
Ucat, _ = concatenate_matrices_with_metadata(AbstractMatrix{AbstractMatrix}(reshape(U_blocks,2,1)))
Vcat, _ = concatenate_matrices_with_metadata(AbstractMatrix{AbstractMatrix}(reshape(V_blocks,2,1)))
U, s, V = svd([A B; C D])
@show Ucat ≈ U
@show Vcat ≈ V
true
true
```

Thus we would like to expect that the U_blocks * S * V_blocks' is approximately equal to the original matrix. The Hermitian conjugate need to be taken into account when dealing with complex matrices. 
This choice of the convension if different from the ITensor package, where the Hermitian conjugate is taken into account in the SVD function.
In VecG tensor svd, we will include the Hermitian conjugate in the SVD function, such that the contraction of U * S * V will directly give the original matrix.
"""
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


function outer_product(A::Array, B::Array)
    return reshape([a * b for a in A, b in B], size(A)..., size(B)...)
end

function test_outer_product()
    A = rand(2,3)
    B = rand(2,2)
    return outer_product(A, B)
end

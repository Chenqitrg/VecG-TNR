using Revise
includet("main.jl")
using .VecG_TNR
using LinearAlgebra

function test_concat_matrix()
    M = Matrix{Matrix}(undef, 2, 2)
    M[1, 1] = [1 2; 3 4]
    M[1, 2] = [5 6; 7 8]
    M[2, 1] = [9 10; 11 12]
    M[2, 2] = [13 14; 15 16]
    @show M
    @show concatenate_matrices_with_metadata(M)

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
    @show metadata
end

function test_svd()
    A = [1 2; 3 4]  # 2x3 matrix
B = [5 6 6; 7 8 8]  # 2x3 matrix
C = [9 10; 11 12; 15 16; 17 18]  # 4x2 matrix
D = [13 14 19; 15 16 20; 21 22 23; 25 26 27]  # 4x3 matrix
matrices = Array{Matrix{Float64}}(undef, 2, 2)
matrices[1, 1] = A
matrices[1, 2] = B
matrices[2, 1] = C
matrices[2, 2] = D
N = [A B; C D]
U_blocks, S, V_blocks = block_matrix_svd(matrices)
U, s, V = svd(N)
display(U)
display(U_blocks[1])
display(U_blocks[2])
display(V)
display(V_blocks[1])
display(V_blocks[2])
display(s)
display(S)
end

function test_svd_complex()
    A = [1+im 2-im; 3+2im 4+3im]  # 2x3 matrix
    B = [5+im 6-im 6+im; 7+2im 8+3im 8+2im]  # 2x3 matrix
    C = [9+im 10-im; 11+2im 12+3im; 15+im 16-im; 17+2im 18+3im]  # 4x2 matrix
    D = [13+im 14-im 19+im; 15+2im 16+3im 20+2im; 21+im 22-im 23+im; 25+2im 26+3im 27+2im]  # 4x3 matrix
    matrices = Array{Matrix{ComplexF64}}(undef, 2, 2)
    matrices[1, 1] = A
    matrices[1, 2] = B
    matrices[2, 1] = C
    matrices[2, 2] = D
    U_blocks, S, V_blocks = block_matrix_svd(matrices)
    N = [A B; C D]
    U, s, V = svd(N)
    @show N' == adjoint(N) # Adjoint can be simplified by a prime
    @show U * diagm(s) * V' ≈ N # ≈ is a symbol for isapprox
    display(U)   
    display(U_blocks[1])  
    display(U_blocks[2])
    display(V)
    display(V_blocks[1])
    display(V_blocks[2])
    display(s)
    display(S)
    @show V[1:2,:] == V_blocks[1]  # true
    @show V[3:5,:] == V_blocks[2]  # true
    @show U[1:2,:] == U_blocks[1]  # true
    @show U[3:6,:] == U_blocks[2]  # true
    @show concatenate_matrices_with_metadata(AbstractMatrix{AbstractMatrix}(reshape(U_blocks,2,1)))[1] == U # true
    @show concatenate_matrices_with_metadata(AbstractMatrix{AbstractMatrix}(reshape(V_blocks,2,1)))[1] == V # true
end

# test_concat_matrix()
# test_svd()
test_svd_complex()
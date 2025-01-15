include("package/groups.jl")
include("package/block_matrix_calculation.jl")
include("package/VecGtensor.jl")
include("package/VecGdecomposition.jl")
include("package/VecGcontraction.jl")
include("initialization.jl")
include("package/display.jl")

# --1--L--2--leg--mor---- = --leg--mor--1--newL--2--
function single_step_QR(L::Mor{G, T}, mor::Mor{G, T}, leg::Int)  where {T, G<:Group}
    new_mor = VecG_tensordot(L, mor, (2,), (leg,)) # Please be careful that the tensordot will rotate the tensor
    _, newL = VecG_qr(new_mor, (1,2,3))
    return newL
end

#----mor--leg--2--R--1-- = --2--newR--1--mor--leg--
function single_step_LQ(mor::Mor{G,T}, R::Mor{G, T}, leg::Int)  where {T,G<:Group}
    new_mor = VecG_tensordot(mor, R, (leg,), (2,))
    _, newR = VecG_qr(new_mor, (2,3,4))
    return newR
end

#--1--L--2--2--R--1-- = ----1--U--2----1--S--2----2--V--1--
#----2--R--1----2--V†--1----1--S^{-1/2}--2---- = ----1--PR--2----
#----1--S^{-1/2}--2----1--U†--2----1--L--2---- = ----1--PL--2----
function to_projector(L::Mor{G, T}, R::Mor{G, T}, epsilon::Float64) where {T, G<:Group}
    LR = VecG_tensordot(L, R, (2,), (2,))
    U, S, V = VecG_svd(LR, (1,), epsilon)

    PR = VecG_tensordot(VecG_tensordot(R, VecG_dag(V), (1,), (2,)), 1 ./ sqrt.(S), (2,), (1,))
    PL = VecG_tensordot(VecG_tensordot(1 ./sqrt.(S), VecG_dag(U), (2,), (1,)), L, (2,), (1,))
    return PR, PL
end

#         |
#         2
#         |
# ----1---T---3----
#         |
#         4
#         |


function entanglement_filtering(morA::Mor{G, T}, morB::Mor{G, T}, N_ef::Int, epsilon::Float64) where {T, G<:Group}
    obL1 = morA[1]
    obL2 = morB[2]
    obL3 = morA[3]
    obL4 = morB[4]

    obR4 = morB[3]
    obR3 = morA[2]
    obR2 = morB[1]
    obR1 = morA[4]

    L1 = identity_mor(T, obL1)
    L2 = identity_mor(T, obL2)
    L3 = identity_mor(T, obL3)
    L4 = identity_mor(T, obL4)

    R4 = identity_mor(T, obR4)
    R3 = identity_mor(T, obR3)
    R2 = identity_mor(T, obR2)
    R1 = identity_mor(T, obR1)

    for n = 1:N_ef
        L2 = single_step_QR(L1, morA, 1)
        L3 = single_step_QR(L2, morB, 2)
        L4 = single_step_QR(L3, morA, 3)
        L1 = single_step_QR(L4, morB, 4)
        Omega = max_abs(L1)
        L1 = L1 ./ Omega
    end

    for n = 1:N_ef
        R3 = single_step_LQ(morB, R4, 3)
        R2 = single_step_LQ(morA, R3, 2)
        R1 = single_step_LQ(morB, R2, 1)
        R4 = single_step_LQ(morA, R1, 4)
        Omega = max_abs(R4)
        R4 = R4 ./ Omega
    end
    
    P4R, P1L = to_projector(L1, R4, epsilon)
    P1R, P2L = to_projector(L2, R1, epsilon)
    P2R, P3L = to_projector(L3, R2, epsilon)
    P3R, P4L = to_projector(L4, R3, epsilon)

    newmorA = VecG_tensordot(P3R, morA, (1,), (2,)) # Be careful that the legs of morA are rotated whenever contraction is implemented
    newmorA = VecG_tensordot(P3L, newmorA, (2,), (2,))
    newmorA = VecG_tensordot(P1R, newmorA, (1,), (2,))
    newmorA = VecG_tensordot(P1L, newmorA, (2,), (2,)) # We rotate it in four times

    newmorB = VecG_tensordot(P2L, morB, (2,), (2,))
    newmorB = VecG_tensordot(P4R, newmorB, (1,), (2,))
    newmorB = VecG_tensordot(P4L, newmorB, (2,), (2,))
    newmorB = VecG_tensordot(P2R, newmorB, (1,), (2,))

    return newmorA, newmorB
end

function coarse_graining(morA::Mor{G, T}, morB::Mor{G, T}, Dcut::Int, epsilon::Float64) where {T, G<:Group}
    ALU, ARD = VecG_factorize(morA, (1,2), Dcut, epsilon)
    BRU, BLD = VecG_factorize(morB, (2,3), Dcut, epsilon)
    newLU = VecG_tensordot(BRU, ARD, (1,), (2,))
    newRD = VecG_tensordot(BLD, ALU, (1,), (2,))
    newmor = VecG_tensordot(newLU, newRD, (4,1), (4,1))
    Omega = max_abs(newmor)
    return newmor./Omega
end

# function Gilt(morA::Mor{G, T}, morB::Mor{G, T}, epsilon::Float64) where {T, G<:Group}
#     Env = VecG_tensordot(morA, morB, (4,), (2,))
#     Env = VecG_tensordot(Env, morA, (6,), (3,))
#     Env = VecG_tensordot(Env, morB, (8,), (4,))
#     U, S, _ = VecG_svd(Env, (10, 1))
#     t = []

# end

# beta = log(sqrt(2)+1)/2

# mor = ising(1.01 * beta)

# for n = 1:100
#     global mor
#     println("n = $n")
#     morA, morB = entanglement_filtering(mor, mor, 10, 1e-12)
#     mor = coarse_graining(morA, morB, 32, 1e-15)
# end

# @show mor
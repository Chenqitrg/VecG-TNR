include("package/groups.jl")
include("package/VecGtensor.jl")
include("package/VecGdecomposition.jl")
include("package/VecGcontraction.jl")
include("initialization.jl")

function single_step_QR(L::Mor{G, T}, mor::Mor{G, T}, leg::Int)  where {T, G<:Group}
    new_mor = VecG_tensordot(L, mor, (2,), (n,))
    qr_leg = map(x -> mod(x - 1, 4) + 1, (leg, leg+1, leg+2))
    _, newL = VecG_qr(new_mor, qr_leg)
    return newL
end

function single_step_LQ(mor::Mor{G,T}, R::Mor{G, T}, leg::Int)  where {T,G<:Group}
    new_mor = VecG_tensordot(mor, R, (n,), (2,))
    lq_leg = map(x -> mod(x - 1, 4) + 1, (leg-2, leg-1, leg))
    _, newR = VecG_qr(new_mor, lq_leg)
    return newR
end

#       |
#       2
#       |
# --1---T---3--
#       |
#       4
#       |

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
    
    return newA, newB
end
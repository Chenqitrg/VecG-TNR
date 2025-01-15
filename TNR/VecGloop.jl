include("../package/main.jl")
using .VecG_TNR

function to_T_array(TA::Mor{G, T}, TB::Mor{G, T}) where {T, G<:Group}
    return [TA, VecG_permutedims(TB, (2,3,4,1)), VecG_permutedims(TA, (3,4,1,2)), VecG_permutedims(TB, (4,1,2,3))]
end

function to_S_array(loop_T_array::Vector{Mor{G, T}}, D_cut::Int, epsilon::Float64)  where {T, G<:Group}
    loop_S_array = Vector{Mor{G, T}}(undef, 8)

    for site_T = 1 : 4
        loop_S_array[2*site_T - 1], loop_S_array[2*site_T] = VecG_factorize(loop_T_array[site_T], (1,2), D_cut, epsilon)
        loop_S_array[2*site_T] = VecG_permutedims(loop_S_array[2*site_T], (3,1,2))
    end

    return loop_S_array
end

#--2--S†-3--
#     |      
#--1--S--4--

function to_SS_array(loop_S_array::Vector{Mor{G, T}}) where {T, G<:Group}
    loop_SS_array = Vector{Mor{G, T}}(undef, 8)
    for i = 1 : 8
        loop_SS_array[i] = VecG_tensordot(loop_S_array[i], VecG_dag(loop_S_array[i]), (2,), (2,))
#--3--S†-4--
#     |      
#--2--S--1--
    loop_SS_array[i] = VecG_permutedims(loop_SS_array[i], (2,3,4,1))
    end
    return loop_SS_array
end

function to_TT_array(loop_T_array::Vector{Mor{G, T}}) where {T, G<:Group}
    loop_TT_array = Vector{Mor{G, T}}(undef, 4)
    for i = 1 : 4
        loop_TT_array[i] = VecG_tensordot(loop_T_array[i], VecG_dag(loop_T_array[i]), (2,3), (3,2))
#--3--S†-4--
#     |      
#--2--S--1--
        loop_TT_array[i] = VecG_permutedims(loop_TT_array[i], (2,3,4,1))
    end
    return loop_TT_array
end
#---2-S†--S†-3--
#     |   |
#      \ /     
#----1--T--4----

#     |    |
#     2    3
#     |    |
#--1--S----S--4--

function to_TSS_array(loop_T_array::Vector{Mor{G, T}}, loop_S_array::Vector{Mor{G, T}}) where {T, G<:Group}
    loop_TSS_array = Vector{Mor{G, T}}(undef, 4)
    for i = 1 : 4
        SS = VecG_dag(VecG_tensordot(loop_S_array[2*i-1], loop_S_array[2*i], (3,), (1,)))
        loop_TSS_array[i] = VecG_tensordot(loop_T_array[i], SS, (2,3), (3,2))
#---3-S†--S†-4--
#     |   |
#      \ /     
#----2--T--1----
        loop_TSS_array[i] = VecG_permutedims(loop_TSS_array[i], (2,3,4,1))
    end
    return loop_TSS_array
end

function to_number(loop_array::Vector{Mor{G, T}}) where {T, G<:Group}
    cont = loop_array[1]
    for i in 2:length(loop_array)
        cont = VecG_tensordot(cont, loop_array[i], (3,4), (2,1))
    end

    cont = VecG_partial_trace(cont, 2)
    return cont
end

function cost_function(loop_TT_array::Vector{Mor{G, T}}, loop_TSS_array::Vector{Mor{G, T}}, loop_SS_array::Vector{Mor{G, T}}) where {T, G<:Group}
    numTT = to_number(loop_TT_array)
    numTSS = to_number(loop_TSS_array)
    numSS = to_number(loop_SS_array)
    return numTT + numSS - numTSS - conj(numTSS)
end

function loop_initialization(TA::Mor{G, T}, TB::Mor{G, T}, D_cut::Int, epsilon::Float64) where {T, G<:Group}
    loop_T = to_T_array(TA, TB)
    loop_S = to_S_array(loop_T, D_cut, epsilon)
    loop_TT = to_TT_array(loop_T)
    loop_SS = to_SS_array(loop_S)
    loop_TSS = to_TSS_array(loop_T, loop_S)
    return loop_T, loop_S, loop_TT, loop_SS, loop_TSS
end

function to_N(loop_S_array::Vector{Mor{G, T}}, loop_SS_array::Vector{Mor{G, T}}, site::Int) where {T, G<:Group}
    
end

# function to_W(loop_T_array::Vector{ITensor}, loop_S_array::Vector{ITensor}, loop_TSS_array::Vector{ITensor}, site_S::Int64)
#     site_TSS = to_site_T(site_S)
#     run_TSS = to_next_T(site_TSS)

#     W = loop_TSS_array[run_TSS]

#     run_TSS = to_next_T(run_TSS)
#     while run_TSS != site_TSS
#         W *= loop_TSS_array[run_TSS]
#         run_TSS = to_next_T(run_TSS)
#     end

#     T = loop_T_array[site_TSS]
#     ul, ur = get_taggedindices(T, ("up, left", "up, right"))

#     if mod(site_S, 2) == 0
#         S = loop_S_array[site_S - 1]
#         TS = S' * δ(ul', ul) * T * δ(ur, ur')
#         W = TS * W
#     else
#         S = loop_S_array[site_S + 1]
#         TS = S' * δ(ur', ur) * T * δ(ul, ul')
#         W = TS * W
#     end

#     return W, TS
# end

# function to_N(loop_S_array::Vector{ITensor}, loop_SS_array::Vector{ITensor}, site_S::Int64)
#     run_S = to_next_S(site_S)
#     N = loop_SS_array[run_S]

#     run_S = to_next_S(run_S)

#     while run_S != site_S
#         N *= loop_SS_array[run_S]
#         run_S = to_next_S(run_S)
#     end
#     # println("S_contract: $(N * loop_SS_array[site_S])")

#     S = loop_S_array[site_S]
#     u = get_taggedindex(S, "up")

#     N *= δ(u', u)

#     return N
# end

# function update_S(N::ITensor, W::ITensor)
#     a, b, c = inds(W)
#     C_im = combiner(a, b, c, tags = "toim")
#     im = combinedind(C_im)
#     W_1 = W * C_im
#     e, f, g = noprime.((a, b, c))
#     C_dom = combiner(e, f, g, tags = "todom")
#     dom = combinedind(C_dom)
#     N_2 = N * C_im * C_dom

#     W_vec = vector(W_1, im)
#     N_mat = matrix(N_2, im, dom)
#     println("N matrix:")
#     display(N_mat)
#     println("W matrix")
#     display(W_vec)
#     prob = LinearProblem(N_mat, W_vec)
#     sol = solve(prob)
#     Snew_vec = sol.u
#     display(Snew_vec)
#     Snew = ITensor(Snew_vec, dom)
#     Snew *= C_dom
#     println("new S:")
#     # @show Snew
#     return Snew
# end

# function loop_optimization(TA::ITensor, TB::ITensor, D_cut::Int64, N_loop::Int64)
#     TT, loop_T_array, loop_S_array, loop_TSS_array, loop_SS_array = loop_initialization(TA, TB, D_cut)
#     # println("TT = $TT")
#     cost = cost_function(TT, loop_SS_array, loop_TSS_array)
#     println("Cost = $cost")
#     for n = 1 : N_loop
#         println("n_loop = $n \n")
#         for site_S = 1 : 2
#             site_TSS = to_site_T(site_S)
#             # println("original S")
#             # @show loop_S_array[site_S]
#             N = to_N(loop_S_array, loop_SS_array, site_S)
#             W, TS = to_W(loop_T_array, loop_S_array, loop_TSS_array, site_S)
#             S_new = update_S(N, W)
#             u = get_taggedindex(S_new, "up")

#             loop_S_array[site_S] = S_new
#             loop_SS_array[site_S] = S_new' * δ(u', u) * S_new
#             loop_TSS_array[site_TSS] = S_new' * TS 
#             cost = cost_function(TT, loop_SS_array, loop_TSS_array)
#             println("Cost = $cost")
#         end
#     end


#     return loop_S_array
# end
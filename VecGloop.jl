function loop_initialization(TA::ITensor, TB::ITensor, D_cut::Int64)
    loop_T_array = to_loop_T_array(TA, TB)

    loop_S_array = Vector{ITensor}(undef, 8)

    T = loop_T_array[1]
    ul, ur = get_taggedindices(T, ("up, left", "up, right"))
    TT = T' * δ(ul', ul) * δ(ur', ur) * T
    for n = 2 : 4
        T = loop_T_array[n]
        ul, ur = get_taggedindices(T, ("up, left", "up, right"))
        TT_new = T' * δ(ul', ul) * δ(ur', ur) * T
        TT = TT * TT_new
    end

    cost_const = scalar(TT)

    for site_T = 1 : 4
        last_T = to_previous_T(site_T)
        T = loop_T_array[site_T]
        l, ul = get_taggedindices(T, ("link, $last_T, $site_T", "up, left", "up, right"))
        loop_S_array[2*site_T - 1], loop_S_array[2*site_T] = factorize(T, l, ul; which_decomp = "svd", maxdim = D_cut, tags = "link, $site_T")
    end

    loop_TSS_array = Vector{ITensor}(undef, 4)

    for site_T in 1 : 4
        T = loop_T_array[site_T]
        SL, SR = loop_S_array[2 * site_T - 1], loop_S_array[2 * site_T]
        ul, ur = get_taggedindices(T, ("up, left", "up, right"))
        T *= δ(ul', ul) * δ(ur', ur)
        loop_TSS_array[site_T] = SL' * SR' * T 
    end

    loop_SS_array = Vector{ITensor}(undef, 8)

    for site_S in 1 : 8
        S = loop_S_array[site_S]
        u = get_taggedindex(S, "up")
        loop_SS_array[site_S] = S' * δ(u', u) * S
    end

    # println("pass_init")
    return cost_const, loop_T_array, loop_S_array, loop_TSS_array, loop_SS_array
end

function cost_function(cost_const::Float64, loop_SS_array::Vector{ITensor}, loop_TSS_array::Vector{ITensor})
    SS = loop_SS_array[1]
    for n = 2 : 8
        SS *= loop_SS_array[n]
    end
    TSS = loop_TSS_array[1]
    for n = 2 : 4
        TSS *= loop_TSS_array[n]
    end
    cost = cost_const + scalar(SS) - 2 * scalar(TSS)

    return cost
end


function to_W(loop_T_array::Vector{ITensor}, loop_S_array::Vector{ITensor}, loop_TSS_array::Vector{ITensor}, site_S::Int64)
    site_TSS = to_site_T(site_S)
    run_TSS = to_next_T(site_TSS)

    W = loop_TSS_array[run_TSS]

    run_TSS = to_next_T(run_TSS)
    while run_TSS != site_TSS
        W *= loop_TSS_array[run_TSS]
        run_TSS = to_next_T(run_TSS)
    end

    T = loop_T_array[site_TSS]
    ul, ur = get_taggedindices(T, ("up, left", "up, right"))

    if mod(site_S, 2) == 0
        S = loop_S_array[site_S - 1]
        TS = S' * δ(ul', ul) * T * δ(ur, ur')
        W = TS * W
    else
        S = loop_S_array[site_S + 1]
        TS = S' * δ(ur', ur) * T * δ(ul, ul')
        W = TS * W
    end

    return W, TS
end

function to_N(loop_S_array::Vector{ITensor}, loop_SS_array::Vector{ITensor}, site_S::Int64)
    run_S = to_next_S(site_S)
    N = loop_SS_array[run_S]

    run_S = to_next_S(run_S)

    while run_S != site_S
        N *= loop_SS_array[run_S]
        run_S = to_next_S(run_S)
    end
    # println("S_contract: $(N * loop_SS_array[site_S])")

    S = loop_S_array[site_S]
    u = get_taggedindex(S, "up")

    N *= δ(u', u)

    return N
end

function update_S(N::ITensor, W::ITensor)
    a, b, c = inds(W)
    C_im = combiner(a, b, c, tags = "toim")
    im = combinedind(C_im)
    W_1 = W * C_im
    e, f, g = noprime.((a, b, c))
    C_dom = combiner(e, f, g, tags = "todom")
    dom = combinedind(C_dom)
    N_2 = N * C_im * C_dom

    W_vec = vector(W_1, im)
    N_mat = matrix(N_2, im, dom)
    println("N matrix:")
    display(N_mat)
    println("W matrix")
    display(W_vec)
    prob = LinearProblem(N_mat, W_vec)
    sol = solve(prob)
    Snew_vec = sol.u
    display(Snew_vec)
    Snew = ITensor(Snew_vec, dom)
    Snew *= C_dom
    println("new S:")
    # @show Snew
    return Snew
end

function loop_optimization(TA::ITensor, TB::ITensor, D_cut::Int64, N_loop::Int64)
    TT, loop_T_array, loop_S_array, loop_TSS_array, loop_SS_array = loop_initialization(TA, TB, D_cut)
    # println("TT = $TT")
    cost = cost_function(TT, loop_SS_array, loop_TSS_array)
    println("Cost = $cost")
    for n = 1 : N_loop
        println("n_loop = $n \n")
        for site_S = 1 : 2
            site_TSS = to_site_T(site_S)
            # println("original S")
            # @show loop_S_array[site_S]
            N = to_N(loop_S_array, loop_SS_array, site_S)
            W, TS = to_W(loop_T_array, loop_S_array, loop_TSS_array, site_S)
            S_new = update_S(N, W)
            u = get_taggedindex(S_new, "up")

            loop_S_array[site_S] = S_new
            loop_SS_array[site_S] = S_new' * δ(u', u) * S_new
            loop_TSS_array[site_TSS] = S_new' * TS 
            cost = cost_function(TT, loop_SS_array, loop_TSS_array)
            println("Cost = $cost")
        end
    end


    return loop_S_array
end
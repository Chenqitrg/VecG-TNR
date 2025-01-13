using ITensors, LinearSolve, LinearAlgebra

global Tc = 2 / log(sqrt(2) + 1)
global N_filt = 10
global N_loop = 10
global N = 10
global D_cut = 16

function tensor_cZ2(Temporature::Float64)
    β = 1 / Temporature

    s = Index(2)
    l, u, r, d = addtags.(s, ("left", "up", "right", "down"))
    T = ITensor(l, r, u, d)

    Sig(x::Int) = 1 - 2 * (x - 1)

    for sl in 1:2
        for sd in 1:2
            for sr in 1:2
                for su in 1:2
                    E = - (Sig(sl) * Sig(sd) + Sig(sd) * Sig(sr) + Sig(sr) * Sig(su) + Sig(su) * Sig(sl))
                    P = exp(-β * E)
                    T[l=>sl, r=>sr, u=>su, d=>sd] = P
                end
            end
        end
    end
    
    normT = norm(T)
    T /= normT

    return T
end

to_site_T(site_S::Int64) = div(site_S + 1, 2)
to_previous_T(site_T::Int64) = mod(site_T - 2, 4) + 1
to_next_T(site_T::Int64) = mod(site_T, 4) + 1
to_previous_S(site_S::Int64) = mod(site_S - 2, 8) + 1
to_next_S(site_S::Int64) = mod(site_S, 8) + 1

# Given the tag, find the corresponded index 
get_taggedindex(T::ITensor, t::String) = getfirst(x->hastags(x, t), inds(T))
get_taggedindices(T::ITensor, ts::Tuple{Vararg{String}}) = map(x -> get_taggedindex(T, x), ts)
get_taggedindices(T::ITensor) = get_taggedindices(T, ("left", "up", "right", "down"))

function to_projector(L::ITensor, R::ITensor)
    Ll = get_taggedindex(L, "L, left")

    link = commonind(L, R)
    U, S, V = svd(L * R, Ll)

    sqrt_S = sqrt.(S)
    inv_sqrt_S = 1 ./ sqrt_S

    PR = replacetags(R * V * inv_sqrt_S, tags(link)=>"PR, left", "Link,u"=>"PR, right")
    PL = replacetags(inv_sqrt_S * U * L, tags(link)=>"PL, right", "Link,v"=>"PL, left")

    return PR, PL
end

function to_loop_T_array(TA::ITensor, TB::ITensor)
    Al, Au, Ar, Ad = get_taggedindices(TA)
    Bl, Bu, Br, Bd = get_taggedindices(TB)

    link_12, link_23, link_34, link_41, ul1, ur1, ul2, ur2, ul3, ur3, ul4, ur4 = settags.([Bu, Ar, Bd, Al, Au, Ar, Br, Bd, Ad, Al, Bl, Bu], ["link, 1, 2", "link, 2, 3", "link, 3, 4", "link, 4, 1", "up, left, 1", "up, right, 1", "up, left, 2", "up, right, 2", "up, left, 3", "up, right, 3", "up, left, 4", "up, right, 4"])

    loop_T_array = [
        TA * δ(Al, link_41) * δ(Ad, link_12) * δ(Au, ul1) * δ(Ar, ur1);
        TB * δ(Bu, link_12) * δ(Bl, link_23) * δ(Br, ul2) * δ(Bd, ur2);
        TA * δ(Ar, link_23) * δ(Au, link_34) * δ(Ad, ul3) * δ(Al, ur3);
        TB * δ(Bd, link_34) * δ(Br, link_41) * δ(Bl, ul4) * δ(Bu, ur4)
    ]

    return loop_T_array    
end


function entanglement_filtering(TA::ITensor, TB::ITensor, N_filt::Int64)
    Al, Au, Ar, Ad = get_taggedindices(TA)
    Bl, Bu, Br, Bd = get_taggedindices(TB)

    loop_T_array = to_loop_T_array(TA, TB)
    # println("pass T_array")
    # @show inds.(loop_T_array)
    L_array = Vector{ITensor}(undef, 4)
    R_array = Vector{ITensor}(undef, 4)

    Lr = commonind(loop_T_array[4], loop_T_array[1])
    Ll = settags(Lr, "L, left")
    Rl = commonind(loop_T_array[4], loop_T_array[1])
    Rr = settags(Rl, "R, right")
    L = δ(Ll, Lr)
    # println(Lr)
    # println(inds(L))
    L_array[1] = L
    R = δ(Rl, Rr)
    R_array[4] = R

    # println("pass filt setup")

    for n = 1 : N_filt
        # println("n = $n")
        for site_T = 1 : 4
            # println("site_T = $site_T")
            T = loop_T_array[site_T]
            nextsite_T = to_next_T(site_T)
            LT = L * T
            _, L = qr(LT, uniqueinds(LT, loop_T_array[nextsite_T]); tags = "L, left")
            L_array[nextsite_T] = L
        end
        Omega_L = tr(array(L))
        L /= Omega_L
    end

    # println("pass L qr")

    for n = 1 : N_filt
        # println(" = $n")
        for site_T = 4:-1:1
            # println("site_T = $site_T")
            T = loop_T_array[site_T]
            previoussite_T = to_previous_T(site_T)
            RT = R * T
            _, R = qr(RT, uniqueinds(RT, loop_T_array[previoussite_T]); tags = "R, right")
            # println(inds(R))
            R_array[previoussite_T] = R
        end
        Omega_R = tr(array(R))
        R /= Omega_R
    end

    # println("pass R qr")

    PR_array = Vector{ITensor}(undef, 4)
    PL_array = Vector{ITensor}(undef, 4)
    for site_T = 1 : 4
        previous_T = to_previous_T(site_T)
        PR_array[previous_T], PL_array[site_T] = to_projector(L_array[site_T], R_array[previous_T])
    end

    # println("pass projector")
    PR_array_indices = map(x->get_taggedindices(x, ("left", "right")), PR_array)
    PL_array_indices = map(x->get_taggedindices(x, ("left", "right")), PL_array)

    PL_array[1] *= δ(PL_array_indices[1][2], Al)
    PR_array[1] *= δ(PR_array_indices[1][1], Ad)
    PL_array[3] *= δ(PL_array_indices[3][2], Ar)
    PR_array[3] *= δ(PR_array_indices[3][1], Au)
    TA = TA * PL_array[1] * PR_array[1] * PR_array[3] * PL_array[3]

    # println("pass TA")

    PL_array[2] *= δ(PL_array_indices[2][2], Bu)
    PR_array[2] *= δ(PR_array_indices[2][1], Bl)
    PL_array[4] *= δ(PL_array_indices[4][2], Bd)
    PR_array[4] *= δ(PR_array_indices[4][1], Br)
    TB = TB * PL_array[2] * PR_array[2] * PR_array[4] * PL_array[4]

    # println("pass TB")

    Al, Au, Ar, Ad, Bl, Bu, Br, Bd = settags.([PL_array_indices[1][1], PR_array_indices[3][2], PL_array_indices[3][1], PR_array_indices[1][2], PR_array_indices[2][2], PL_array_indices[2][1], PR_array_indices[4][2], PL_array_indices[4][1]], ["A, left", "A, up", "A, right", "A, down", "B, left", "B, up", "B, right", "B, down"])

    TA = TA * δ(PL_array_indices[1][1], Al) * δ(PL_array_indices[3][1], Ar) * δ(PR_array_indices[3][2], Au) * δ(PR_array_indices[1][2], Ad)
    TB = TB * δ(PR_array_indices[2][2], Bl) * δ(PR_array_indices[4][2], Br) * δ(PL_array_indices[2][1], Bu) * δ(PL_array_indices[4][1], Bd)

    # println("pass contraction")
    
    return TA, TB
end

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

function contract_loop(loop_S_array::Vector{ITensor})
    Aru, Bru, Brd, Ard, Ald, Bld, Blu, Alu = get_taggedindex.(loop_S_array, "up")

    TA = loop_S_array[8] * loop_S_array[1] * δ(Aru, Ard) * δ(Alu, Ald) * loop_S_array[4] * loop_S_array[5]
    TB = loop_S_array[2] * loop_S_array[3] * δ(Bru, Blu) * δ(Brd, Bld) * loop_S_array[6] * loop_S_array[7]

    TB = replacetags!(TB, "link, 2"=>"B, left", "link, 1"=>"B, up", "link, 4"=>"B, right", "link, 3"=>"B, down")
    TA = replacetags!(TA, "link, 4"=>"A, left", "link, 3"=>"A, up", "link, 2"=>"A, right", "link, 1"=>"A, down")

    Al, Au, Ar, Ad = get_taggedindices(TA)
    Z = scalar(TA * δ(Al, Ar) * δ(Au, Ad))
    println("Z = $Z \n")
    TA /= Z
    TB /= Z

    return TA, TB
end

function conformal_spectrum(TA::ITensor, TB::ITensor)
    Al, Au, Ar, Ad = get_taggedindices(TA)
    Bl, Bu, Br, Bd = get_taggedindices(TB)
    P12 = TA * δ(Ad, Bu) * TB * δ(Bd, Au)
    U = combiner(Al, Bl)
    D = combiner(Br', Ar')
    P = matrix(P12 * δ(Ar, Bl') * δ(Br, Al') * P12' * U * D)
    # println("Proj")
    spec22_v = sort(abs.(eigvals(P)), rev = true)
    # println("ev")
    h = - log.(spec22_v ./ spec22_v[1]) / (2*pi)

    return h
end

T = tensor_cZ2(Tc)

TA, TB = addtags.([T, T], ["A", "B"])

for n = 1 : 1
    global TA, TB

    println("n = $n")
    TA, TB = entanglement_filtering(TA, TB, 5)
    # println("pass_filt")

    # local loop_S_array = loop_optimization(TA, TB, 16, 1)

    loop_S_array = loop_optimization(TA, TB, 16, 1)
    # TT, loop_T_array, loop_S_array, loop_TSS_array, loop_SS_array = loop_initialization(TA, TB, D_cut)

    # for nloop = 1 : 5
    #     println("n_loop = $nloop")
    #     for site_S = 1 : 8
    #         site_TSS = to_site_T(site_S)

    #         N = to_N(loop_S_array, loop_SS_array, site_S)
    #         W, TS = to_W(loop_T_array, loop_S_array, loop_TSS_array, site_S)
    #         S_new = update_S(N, W)
    #         u = get_taggedindex(S_new, "up")

    #         loop_S_array[site_S] = S_new 
    #         loop_SS_array[site_S] = S_new' * δ(u', u) * S_new
    #         loop_TSS_array[site_TSS] = S_new' * TS 
    #     end
    # end

    TA, TB = contract_loop(loop_S_array)
end

# display(TA)
conformal_spectrum(TA, TB)
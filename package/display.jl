# Custom display for a group
function Base.show(io::IO, G::Group)
    if G isa CyclicGroup
        print(io, "ℤ", subscript(G.n))
    elseif G isa DihedralGroup
        print(io, "𝔻", subscript(2 * G.n))
    end
end

# Custom display for a group element in a CyclicGroup
function Base.show(io::IO, x::GroupElement{CyclicGroup})
    if x.value == 0
        print(io, "e")
    elseif x.value == 1
        print(io, "a")
    else
        print(io, "a", superscript(x.value))
    end
end

# Custom display for a group element in a DihedralGroup
function Base.show(io::IO, x::GroupElement{DihedralGroup})
    s, r = x.value
    if s == 0 && r == 0
        print(io, "e")
    elseif s == 0 && r == 1
        print(io, "r")
    elseif s == 0 && r > 0
        print(io, "r", superscript(r))
    elseif s == 1 && r == 0
        print(io, "s")
    elseif s == 1 && r == 1
        print(io, "sr")
    else
        print(io, "sr", superscript(r))
    end
end


function subscript(n::Integer)
    # 定义一个映射，将字符 '0'-'9' 转换为 Unicode 下角标字符
    subs = Dict('0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄',
                '5' => '₅', '6' => '₆', '7' => '₇', '8' => '₈', '9' => '₉')
    # 将数字 n 转换为字符串，逐位映射到下角标字符，并连接成新的字符串
    return join([subs[c] for c in string(n)])
end

function to_subscript(string::String)
    # 映射表：普通字母 -> 下标字母
    subscript_map = Dict(
        'a' => 'ₐ',
        'e' => 'ₑ',
        'h' => 'ₕ',
        'i' => 'ᵢ',
        'j' => 'ⱼ',
        'k' => 'ₖ',
        'l' => 'ₗ',
        'm' => 'ₘ',
        'n' => 'ₙ',
        'o' => 'ₒ',
        'p' => 'ₚ',
        'r' => 'ᵣ',
        's' => 'ₛ',
        't' => 'ₜ',
        'u' => 'ᵤ',
        'v' => 'ᵥ',
        'x' => 'ₓ'
    )
    # 返回对应的下标字母，如果没有对应的，返回原字母
    return join([subscript_map[c] for c in string])
end

function superscript(n::Integer)
    # 定义一个映射，将字符 '0'-'9' 转换为 Unicode 下角标字符
    subs = Dict('0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
                '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹')
    # 将数字 n 转换为字符串，逐位映射到下角标字符，并连接成新的字符串
    return join([subs[c] for c in string(n)])
end

# Custom display for Object
function Base.show(io::IO, x::Obj)
    group = get_group(x)
    sumd = x.sumd
    el = elements(group)
    op = false
    for i in el
        if sumd[i] != 0
            if op == true
                print(io, " ⊕ ")
            end
            if sumd[i] == 1
                print(io, i)
            else
                print(io, sumd[i], i)
            end
            op = true
        end
    end
end

# Custom display for a sector
function Base.show(io::IO, S::Sector)
    for (i, g) in enumerate(S.sect)
        if i != 1
            print(io, " ⊗ ")
        end
        print(io, g)
    end
end

# Custom display for Morphism
function Base.show(io::IO, ::MIME"text/plain", T::Mor)
    println(io, "Group: ", get_group(T))
    for (i, obj) in enumerate(T.objects)
        println(io, "Leg $i is of object $obj")
    end
end

# Custom display for Morphism
function Base.show(io::IO, ::MIME"text/plain", T::TubeMor)
    println(io, "Group: ", get_group(T))
    for (i, obj) in enumerate(T.in)
        println(io, "In legs: leg $i is of object $obj")
    end
    for (i, obj) in enumerate(T.out)
        println(io, "Out legs: leg $i is of object $obj")
    end
end

# Custom display for a sector in Tube category
function Base.show(io::IO, S::TubeSector)
    for (i, g) in enumerate(S.insect)
        if i != 1
            print(io, " ⊗ ")
        end
        print(io, g)
    end
    print(io, " →(", S.circ, ") ")
    # print(io, sect[2])
    # print(io, "\ ")
    for (i, g) in enumerate(S.outsect)
        if i != 1
            print(io, " ⊗ ")
        end
        print(io, g)
    end
end
struct TubeSector{G<:Group}
    insect::Tuple{Vararg{GroupElement{G}}}
    circ::GroupElement{G}
    outsect::Tuple{Vararg{GroupElement{G}}}
    function TubeSector(in_tup::Tuple{Vararg{GroupElement{G}}}, circ::GroupElement{G}, out_tup::Tuple{Vararg{GroupElement{G}}}) where G<:Group
        in_g = multiply(in_tup)
        out_g = multiply(out_tup)
        if circ * in_g == out_g * circ
            return new{G}(in_tup, circ, out_tup)
        else
            throw(ArgumentError("The sector $in_tup, $circ, $out_tup is invalid"))
        end
    end
end

struct TubeMor{G<:Group, T}
    in::Tuple{Vararg{Obj{G}}}
    out::Tuple{Vararg{Obj{G}}}
    data::Dict{TubeSector{G}, Array{T}}
end

# Construction function for a morphism
function TubeMor(element_type::Type, in::Tuple{Vararg{Obj{G}}}, out::Tuple{Vararg{Obj{G}}}) where {G<:Group}
    data = Dict{TubeSector{G}, Array{element_type}}()
    return TubeMor{G,element_type}(in, out, data)
end

# Get the tensor of the S sector
function Base.getindex(mor::TubeMor{G, T}, S::Sector{G}) where {T, G<:Group}
    return mor.data[S]
end

# Get the tensor of the S sector
function Base.getindex(mor::TubeMor{G, T}, g...) where {T, G<:Group}
    return mor.data[TubeSector(g[1], g[2], g[3])]
end

# Get the group of a morphism
function get_group(T::TubeMor)
    return get_group(T.in[1])
end

function Base.getindex(mor::TubeMor{G, T}, S::Sector{G}) where {T, G<:Group}
    return mor.data[S]
end

# Get the size of the tensor for a given sector
function get_sector_size(mor::TubeMor{G, T}, sector::TubeSector{G}) where {T, G<:Group}
    in_sect = sector.insect
    out_sect = sector.outsect
    in_size = ()
    for (i, g) in enumerate(in_sect)
        object = mor.in[i]
        multiplicity = object[g]
        in_size = (in_size..., multiplicity)
    end

    out_size = ()
    for (i, g) in enumerate(out_sect)
        object = mor.out[i]
        multiplicity = object[g]
        out_size = (out_size..., multiplicity)
    end
    return in_size, out_size
end

# Set a specific sector
function Base.setindex!(mor::TubeMor{G, T}, value::Array{T}, S::TubeSector{G}) where {T, G <: Group}
    terminal_size_in, terminal_size_out = get_sector_size(mor, S)
    data_size = size(value)
    if (terminal_size_in..., terminal_size_out...) != data_size
        throw(ArgumentError("Data size $data_size does not match sector size $terminal_size"))
    end
    mor.data[S] = value
end

# Construct a random morphism in tube category for a given objects
function random_mor(element_type::Type, in::Tuple{Vararg{Obj}}, out::Tuple{Vararg{Obj}})
    mor = TubeMor(element_type, in, out)
    group = get_group(mor)
    for g_circ in elements(group)
        for groups_in in group_iter(group, length(in))
            out_sect = (g_circ * multiply(groups_in)) * inverse(g_circ)
            for groups_out in group_tree(out_sect, length(out))
                sector = TubeSector(groups_out, g_circ, groups_in)
                sector_size_in, sector_size_out = get_sector_size(mor, sector)
                mor[sector] = rand(element_type, (sector_size_in..., sector_size_out...))
            end
        end
    end

    return mor
end

# Construct a random morphism in tube category for a given objects
function zero_mor(element_type::Type, in::Tuple{Vararg{Obj}}, out::Tuple{Vararg{Obj}})
    mor = TubeMor(element_type, in, out)
    group = get_group(mor)
    for g_circ in elements(group)
        for groups_in in group_iter(group, length(in))
            out_sect = (g_circ * multiply(groups_in)) * inverse(g_circ)
            for groups_out in group_tree(out_sect, length(out))
                sector = TubeSector(groups_out, g_circ, groups_in)
                sector_size_in, sector_size_out = get_sector_size(mor, sector)
                mor[sector] = zeros(element_type, (sector_size_in..., sector_size_out...))
            end
        end
    end

    return mor
end


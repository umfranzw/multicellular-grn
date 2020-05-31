module ProteinStoreMod

using ProteinMod
using ProteinPropsMod

export ProteinStore

mutable struct ProteinStore
    proteins::Dict{ProteinPropsMod.ProteinType, Dict{ProteinProps, Protein}}
    
    function ProteinStore()
        new(build_dict())
    end
end

function build_dict()
    proteins = Dict{ProteinPropsMod.ProteinType, Dict{ProteinProps, Protein}}()
    for type in instances(ProteinPropsMod.ProteinType)
        proteins[type] = Dict{ProteinProps, Protein}()
    end

    proteins
end

function num_proteins(ps::ProteinStore)
    total = 0
    for type in keys(ps.proteins)
        total += length(ps.proteins[type])
    end

    total
end

function clear(ps::ProteinStore)
    ps.proteins = build_dict()
end

function contains(ps::ProteinStore, protein::Protein)
    types = instances(ProteinPropsMode.ProteinType)
    i = 0
    found = false
    while !found && i < length(types)
        found = protein.props in keys(ps.proteins[types[i]])
        i += 1
    end

    found
end

function insert(ps::ProteinStore, protein::Protein)
    sub_dict = ps.proteins[protein.props.type]
    sub_dict[protein.props] = protein
end

function alter_concs(ps::ProteinStore, delta::Array{Float64, 1})
    protein = ps.proteins[protein.props.type][protein.props]
    protein.concs = clamp.(protein.concs .+ delta, 0.0, 1.0)
end

function remove(ps::ProteinStore, protein::Protein)
    sub_dict = ps.proteins[protein.props.type]
    if protein.props in keys(sub_dict)
        delete!(sub_dict, protein.props)
    end
end

function get(ps::ProteinStore, props::ProteinProps)
    result = nothing
    sub_dict = ps.proteins[props.type]
    if props in keys(sub_dict)
        result = sub_dict[props]
    end

    result
end

function get_neighbour_proteins_by_loc(ps::ProteinStore, loc::ProteinPropsMod.ProteinLoc, src_cell_id::UInt64)
    neighbours = Array{Protein, 1}()
    for protein in values(ps.proteins[ProteinPropsMod.Neighbour])
        if protein.props.loc == loc && protein.src_cell_id == src_cell_id
            push!(neighbours, protein)
        end
    end

    neighbours
end

function get_by_type(ps::ProteinStore, type::ProteinPropsMod.ProteinType)
    values(ps.proteins[type])
end

function get_by_types(ps::ProteinStore, types::Set{ProteinPropsMod.ProteinType})
    result = Array{Protein, 1}()
    for type in types
        append!(result, values(ps.proteins[type]))
    end

    result
end

function get_all(ps::ProteinStore)
    proteins = Array{Protein, 1}()
    for type in instances(ProteinPropsMod.ProteinType)
        append!(proteins, values(ps.proteins[type]))
    end

    proteins
end

function get_all_props(ps::ProteinStore)
    props = Set{ProteinProps}()
    for type in instances(ProteinPropsMod.ProteinType)
        union!(props, keys(ps.proteins[type]))
    end
    
    props
end

end

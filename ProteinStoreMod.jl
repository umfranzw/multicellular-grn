module ProteinStoreMod

using ProteinMod
using ProteinPropsMod
using RunMod

export ProteinStore

mutable struct ProteinStore
    config::Config
    proteins::Dict{ProteinPropsMod.ProteinType, Dict{ProteinProps, Protein}}
    #note: these are also stored in proteins (above) - this is just for efficiency since they'll be accessed by location a lot
    loc_to_sensor::Dict{ProteinPropsMod.ProteinLoc, Dict{ProteinProps, Protein}}
    
    function ProteinStore(config::Config)
        new(config, build_dicts()...)
    end
end

function build_dicts()
    proteins = Dict{ProteinPropsMod.ProteinType, Dict{ProteinProps, Protein}}()
    
    for type in instances(ProteinPropsMod.ProteinType)
        proteins[type] = Dict{ProteinProps, Protein}()
    end
    
    loc_to_sensor = Dict{ProteinPropsMod.ProteinLoc, Dict{ProteinProps, Protein}}()
    for loc in instances(ProteinPropsMod.ProteinLoc)
        loc_to_sensor[loc] = Dict{ProteinProps, Protein}()
    end

    (proteins, loc_to_sensor)
end

function clear(ps::ProteinStore)
    ps.proteins, ps.loc_to_sensor = build_dicts()
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

    if protein.props.action == ProteinPropsMod.Sensor
        ps.loc_to_sensor[protein.props.loc][protein.props] = protein
    end
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

    if protein.action == ProteinPropsMod.Sensor
        delete!(ps.loc_to_sensor[protein.props.loc], protein.props)
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

function get_by_type(ps::ProteinStore, type::ProteinPropsMod.ProteinType)
    values(ps.proteins[type])
end

function get_sensors_at(ps::ProteinStore, loc::ProteinPropsMod.ProteinLoc)
    values(ps.loc_to_sensor[loc])
end

function get_all(ps::ProteinStore)
    proteins = Array{Protein, 1}()
    for type in instances(ProteinPropsMod.ProteinType)
        push!(proteins, values(ps.proteins[type])...)
    end

    proteins
end

end

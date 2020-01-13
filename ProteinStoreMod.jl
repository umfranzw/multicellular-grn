module ProteinStoreMod

using ProteinMod
using ProteinPropsMod
using RunMod

export ProteinStore

mutable struct ProteinStore
    config::Config
    proteins::Dict{ProteinPropsMod.ProteinTarget, Dict{ProteinProps, Protein}}
    owned_intercell_proteins::Set{ProteinProps} #set of inter-cell proteins that the cell that owns this store has produced
    
    function ProteinStore(config::Config)
        proteins = Dict{ProteinPropsMod.ProteinTarget, Dict{ProteinProps, Protein}}()
        owned_intercell_proteins = Set{ProteinProps}()
        
        for target in instances(ProteinPropsMod.ProteinTarget)
            proteins[target] = Dict{ProteinProps, Protein}()
        end
        
        new(config, proteins, owned_intercell_proteins)
    end
end

function clear(ps::ProteinStore)
    ps.proteins = Dict{ProteinPropsMod.ProteinTarget, Dict{ProteinProps, Protein}}()
    ps.owned_intercell_proteins = Set{ProteinProps}()
    
    for target in instances(ProteinPropsMod.ProteinTarget)
        ps.proteins[target] = Dict{ProteinProps, Protein}()
    end
end

function contains(ps::ProteinStore, protein::Protein)
    for target in instances(ProteinPropsMod.ProteinTarget)
        if protein.props in keys(ps.proteins[target])
            return true
        end
    end

    return false
end

#owned should be set to true if the cell that's inserting the protein produced it
function insert(ps::ProteinStore, protein::Protein, owned::Bool)
    sub_dict = ps.proteins[protein.props.target]
    if protein.props in keys(sub_dict)
        #add new protein concs to existing ones
        for i in 1:length(protein.concs)
            sub_dict[protein.props].concs[i] = min(sub_dict[protein.props].concs[i] + protein.concs[i], 1.0)
        end
    else
        sub_dict[protein.props] = protein
    end

    if owned && protein.props.target == ProteinPropsMod.Inter
        #note: no need to check if it's already present since this is a Set
        push!(ps.owned_intercell_proteins, protein.props)
    end
end

function remove(ps::ProteinStore, protein::Protein)
    sub_dict = ps.proteins[protein.props.target]
    if protein.props in keys(sub_dict)
        delete!(sub_dict, protein.props)
    end
    
    if protein.props.target == ProteinPropsMod.Inter && protein.props in ps.owned_intercell_proteins
        delete!(ps.owned_intercell_proteins, protein.props)
    end
end

function is_owned_intercell_protein(ps::ProteinStore, protein::Protein)
    protein.props in ps.owned_intercell_proteins
end

function get_owned_intercell_proteins(ps::ProteinStore, protein::Protein)
    map(props -> ps.proteins[props.target][props], ps.owned_intercell_proteins)
end

function get(ps::ProteinStore, props::ProteinProps)
    sub_dict = ps.proteins[props.target]
    if props in keys(sub_dict)
        return sub_dict[props]
    else
        return nothing
    end
end

function get_by_target(ps::ProteinStore, target::ProteinPropsMod.ProteinTarget)
    values(ps.proteins[target])
end

function get_by_type(ps::ProteinStore, type::ProteinPropsMod.ProteinType)
    proteins = Array{Protein, 1}()
    for target in instances(ProteinPropsMod.ProteinTarget)
        for protein in values(ps.proteins[target])
            if protein.props.type == type
                push!(proteins, protein)
            end
        end
    end

    proteins
end

function get_all(ps::ProteinStore)
    proteins = Array{Protein, 1}()
    for target in instances(ProteinPropsMod.ProteinTarget)
        push!(proteins, values(ps.proteins[target])...)
    end

    proteins
end

end

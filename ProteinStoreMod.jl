module ProteinStoreMod

using ProteinMod
using RunMod

export ProteinStore,
    insert_protein, get_protein, insert_cell

mutable struct ProteinStore
    run::Run
    proteins::Dict{ProteinMod.ProteinTarget, Dict{ProteinProps, Protein}}
    owned_intercell_proteins::Set{ProteinProps} #set of inter-cell proteins that the cell that owns this store has produced
    
    function ProteinStore(run::Run)
        proteins = Dict{ProteinMod.Scope, Dict{ProteinProps, Protein}}()
        owned_intercell_proteins = Set{ProteinProps}()
        
        for target in instances(ProteinMod.ProteinTarget)
            proteins[target] = Dict{ProteinProps, Protein}()
        end
        
        new(run, proteins, owned_intercell_proteins)
    end
end

function contains(ps::ProteinStore, protein::Protein)
    i = 1
    found = false
    while !found && i <= length(instances(ProteinMod.ProteinTarget))
        target = ProteinMod.ProteinTarget(i)
        found = protein.props in keys(ps.proteins[target])
    end

    found
end

#owned should be set to true if the cell that's inserting the protein produced it
function insert(ps::ProteinStore, protein::Protein, owned::Bool)
    sub_dict = ps.proteins[protein.props.target]
    if protein.props in keys(sub_dict)
        #add new protein concs to existing ones
        for i in 1:length(protein.concs)
            sub_dict[protein.seq].concs[i] = min(sub_dict[protein.seq].concs[i] + protein.concs[i], 1.0)
        end
    else
        sub_dict[protein.props] = protein
    end

    if owned && ProteinMod.get_target(protein) == ProteinMod.Inter
        #note: no need to check if it's already present since this is a Set
        push!(ps.owned_intercell_proteins, protein.props)
    end
end

function is_owned_intercell_protein(ps::ProteinStore, protein::Protein)
    protein.props in ps.owned_intercell_proteins
end

function get(ps::ProteinStore, props::ProteinProps)
    sub_dict = ps.proteins[props.target]
    if props in keys(sub_dict)
        return sub_dict[props]
    else
        return nothing
    end
end

function get_by_target(ps::ProteinStore, target::ProteinMod.ProteinTarget)
    values(ps.proteins[target])
end

function get_by_type(ps::ProteinStore, type::ProteinMod.ProteinType)
    proteins = []
    for target in instances(ProteinMod.ProteinTarget)
        for protein in values(ps.proteins[target])
            if protein.props.type == type
                append!(proteins, protein)
            end
        end
    end

    proteins
end

end

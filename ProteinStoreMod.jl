module ProteinStoreMod

using ProteinMod
using RunMod

export ProteinStore,
    insert_protein, get_protein, insert_cell

mutable struct ProteinStore
    run::Run
    proteins::Dict{ProteinMod.Scope, Dict{BitArray{1}, Protein}}
    owned_intercell_proteins::Set{BitArray{1}}
    
    function ProteinStore(run::Run)
        proteins = Dict{ProteinMod.Scope, Dict{BitArray{1}, Protein}}()
        owned_intercell_proteins = Set{BitArray{1}}()
        
        for scope in instances(ProteinMod.Scope)
            proteins[scope] = Dict{BitArray{1}, Protein}()
        end
        
        new(run, proteins, owned_intercell_proteins)
    end
end

function insert_protein(ps::ProteinStore, protein::Protein, owned::Bool)
    sub_dict = ps.proteins[ProteinMod.get_scope(protein)]
    if protein.seq in keys(sub_dict)
        #add new protein concs to existing ones
        for i in 1:length(protein.concs)
            sub_dict[protein.seq].concs[i] = min(sub_dict[protein.seq].concs[i] + protein.concs[i], 1.0)
        end
    else
        sub_dict[protein.seq] = protein
    end

    if owned && ProteinMod.get_scope(protein) == ProteinMod.InterCell
        #note: no need to check if it's already present since this is a Set
        push!(ps.owned_intercell_proteins, protein.seq)
    end
end

function is_owned_intercell_protein(ps::ProteinStore, protein::Protein)
    protein.seq in ps.owned_intercell_proteins
end

function get_protein(ps::ProteinStore, seq::BitArray{1})
    sub_dict = ps.proteins[ProteinMod.get_scope(seq)]
    if seq in keys(sub_dict)
        return sub_dict[seq]
    else
        return nothing
    end
end

function get_proteins_by_scope(ps::ProteinStore, scope::ProteinMod.Scope)
    values(ps.proteins[scope])
end

function get_proteins_by_target(ps::ProteinStore, target::ProteinMod.Target)
    proteins = []
    for scope in instances(ProteinMod.Scope)
        for protein in values(ps.proteins[scope])
            if ProteinMod.get_target(protein) == target
                append!(proteins, protein)
            end
        end
    end
end

function get_proteins_by_prop(ps::ProteinStore, scope::ProteinMod.Scope, target::ProteinMod.Target)
    proteins = values(ps.proteins[scope])
    filter(p -> ProteinMod.get_target(p) == target, proteins)
end

end

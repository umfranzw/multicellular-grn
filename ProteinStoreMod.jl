module ProteinStoreMod

using ProteinMod
using RunMod

export ProteinStore,
    insert_protein, get_protein, insert_cell

mutable struct ProteinStore
    run::Run
    proteins::Dict{ProteinMod.Scope, Dict{BitArray{1}, Protein}}

    function ProteinStore(run::Run)
        proteins = Dict{ProteinMod.Scope, Dict{BitArray{1}, Protein}}()
        for scope in instances(ProteinMod.Scope)
            proteins[scope] = Dict{BitArray{1}, Protein}()
        end
        
        new(run, proteins)
    end
end

function insert_protein(ps::ProteinStore, protein::Protein)
    sub_dict = ps.proteins[ProteinMod.get_scope(protein)]
    if protein.seq in keys(sub_dict)
        #add new protein concs to existing ones
        for i in 1:length(protein.concs)
            sub_dict[protein.seq].concs[i] = min(sub_dict[protein.seq].concs[i] + protein.concs[i], 1.0)
        end
    else
        sub_dict[protein.seq] = protein
    end
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
    ps.proteins[scope]
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

function insert_cell(ps::ProteinStore, cell_index::Int64)
    for scope in instances(ProteinMod.Scope)
        for protein in values(ps.proteins[scope])
            insert!(protein.concs, cell_index, zeros(Float64, ps.run.num_genes))
        end
    end
end

end

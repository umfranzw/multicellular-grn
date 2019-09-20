module ProteinStoreMod

using ProteinMod
using RunMod

export ProteinStore,
    insert_protein, get_protein, insert_cell

mutable struct ProteinStore
    run::Run
    proteins::Dict{BitArray{1}, Protein}

    function ProteinStore(run::Run)
        new(run, Dict{BitArray{1}, Protein}())
    end
end

function insert_protein(ps::ProteinStore, protein::Protein)
    if protein.seq in keys(ps.proteins)
        #add new protein concs to existing ones
        for i in 1:length(protein.concs)
            ps[protein.seq].concs[i] = min(ps[protein.seq].concs[i] + protein.concs[i], 1.0)
        end
    else
        ps.proteins[protein.seq] = protein
    end
end

function get_protein(ps::ProteinStore, seq::BitArray{1})
    try
        return ps.proteins[seq]
    catch err
        if isa(err, KeyError)
            return nothing
        else
            throw(err)
        end
    end
end

function insert_cell(ps::ProteinStore, cell_index::Int64)
    for protein in values(ps.proteins)
        insert!(protein.concs, zeros(Float64, ps.run.num_genes))
    end
end

end

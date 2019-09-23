module ProteinStoreMod

using ProteinMod
using CellMod
using RunMod

export ProteinStore,
    insert_protein, get_protein, insert_cell

mutable struct ProteinStore
    run::Run
    proteins::Dict{ProteinMod.Scope, Dict{BitArray{1}, Protein}}
    src_cells::Dict{BitArray{1}, Set{Cell}}
    
    function ProteinStore(run::Run)
        proteins = Dict{ProteinMod.Scope, Dict{BitArray{1}, Protein}}()
        src_cells = Dict{BitArray{1}, Set{Cell}}()
        
        for scope in instances(ProteinMod.Scope)
            proteins[scope] = Dict{BitArray{1}, Protein}()
        end
        
        new(run, proteins, src_cells)
    end
end

function insert_protein(ps::ProteinStore, protein::Protein, src_cell::Union{Cell, Nothing})
    sub_dict = ps.proteins[ProteinMod.get_scope(protein)]
    if protein.seq in keys(sub_dict)
        #add new protein concs to existing ones
        for i in 1:length(protein.concs)
            sub_dict[protein.seq].concs[i] = min(sub_dict[protein.seq].concs[i] + protein.concs[i], 1.0)
        end
    else
        sub_dict[protein.seq] = protein
    end

    if src_cell != nothing
        if protein.seq ∉ keys(ps.src_cells)
            ps.src_cells[protein.seq] = Set{Cell}()
        end
        
        #note: no need to check if it's already present since this is a Set
        push!(ps.src_cells[protein.seq], src_cell)
    end
end

function get_src_cells(ps::ProteinStore, protein::Protein)
    if protein.seq in keys(ps.src_cells)
        return ps.src_cells[protein.seq]
    else
        return Set{Cell}()
    end
end

function has_src_cell(ps::ProteinStore, protein::Protein, src_cell::Cell)
    protein.seq in keys(ps.src_cells) && src_cell in ps.src_cells[protein.seq]
end

function add_src_cell(ps::ProteinStore, protein, src_cell::Cell)
    if protein.seq in keys(ps.src_cells)
        push!(ps.src_cells[protein.seq], src_cell)
    else
        error("Error (ProteinStoreMod): add_src_cell() called with non-existent protein.")
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

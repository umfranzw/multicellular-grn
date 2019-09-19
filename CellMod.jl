module CellMod

using ProteinMod
using RunMod
using GeneStateMod
using GeneMod

export Cell

mutable struct Cell
    run::Run
    proteins::Dict{ProteinScope, Dict{BitArray, Protein}}
    gene_states::Array{GeneState, 1}

    function Cell(run::Run, genes::Array{Gene, 1}, initial_proteins::Array{Protein, 1})
        proteins = Dict{ProteinScope, Dict{BitArray, Protein}}()
        proteins[ProteinMod.IntraCell] = Dict{BitArray, Protein}()
        proteins[ProteinMod.InterCell] = Dict{BitArray, Protein}()
        
        for protein in initial_proteins
            scope = ProteinMod.get_scope(protein)
            proteins[scope][protein.seq] = protein
        end

        gene_states = map(g -> GeneState(run, g), genes)

        new(run, proteins, gene_states)
    end
end

end

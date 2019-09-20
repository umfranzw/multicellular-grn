module CellMod

using ProteinMod
using RunMod
using GeneStateMod
using GeneMod
using ProteinStoreMod

export Cell

mutable struct Cell
    run::Run
    gene_states::Array{GeneState, 1}
    store::ProteinStore

    function Cell(run::Run, genes::Array{Gene, 1}, initial_proteins::Array{Protein, 1}, store::ProteinStore)
        for protein in initial_proteins
            ProteinStoreMod.insert_protein(store, protein)
        end

        gene_states = map(g -> GeneState(run, g), genes)

        new(run, gene_states, store)
    end
end

end

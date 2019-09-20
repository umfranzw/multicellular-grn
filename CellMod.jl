module CellMod

using RunMod
using GeneStateMod
using GeneMod

export Cell

mutable struct Cell
    run::Run
    gene_states::Array{GeneState, 1}

    function Cell(run::Run, genes::Array{Gene, 1})
        gene_states = map(g -> GeneState(run, g), genes)

        new(run, gene_states)
    end
end

end

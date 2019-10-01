module CellMod

using RunMod
using GeneStateMod
using GeneMod
using ProteinStoreMod
using SymMod: Sym

export Cell

mutable struct Cell
    run::Run
    gene_states::Array{GeneState, 1}
    proteins::ProteinStore
    energy::Float64
    parent::Union{Cell, Nothing}
    children::Array{Cell, 1}
    sym::Union{Sym, Nothing}

    function Cell(run::Run, genes::Array{Gene, 1}, parent::Union{Cell, Nothing}, sym::Union{Sym, Nothing})
        gene_states = map(g -> GeneState(run, g), genes)
        proteins = ProteinStore(run) #all proteins present in this cell
        children = Array{Cell, 1}()

        new(run, gene_states, proteins, run.initial_cell_energy, nothing, children, sym)
    end

    function Cell(run::Run, gene_states::Array{GeneState, 1}, parent::Union{Cell, Nothing}, sym::Union{Sym, Nothing})
        new(run, gene_states, ProteinStore(), run.initial_cell_energy, Array{Cell, 1}(), parent, sym)
    end
end

end

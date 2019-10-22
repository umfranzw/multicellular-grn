module CellMod

using RunMod
using GeneStateMod
using GeneMod
using ProteinStoreMod
using SymMod: Sym
using MiscUtilsMod

import Base.show

export Cell

mutable struct Cell
    run::Run
    gene_states::Array{GeneState, 1}
    proteins::ProteinStore
    energy::Float64
    parent::Union{Cell, Nothing}
    children::Array{Cell, 1}
    sym::Union{Sym, Nothing}

    function Cell(run::Run, genes::Array{Gene, 1}, sym::Union{Sym, Nothing})
        gene_states = map(g -> GeneState(run, g), genes)
        proteins = ProteinStore(run) #all proteins present in this cell
        children = Array{Cell, 1}()

        cell = new(run, gene_states, proteins, run.initial_cell_energy, nothing, children, sym)
        
        cell
    end

    function Cell(run::Run, gene_states::Array{GeneState, 1}, sym::Union{Sym, Nothing})
        cell = new(run, gene_states, ProteinStore(), run.initial_cell_energy, nothing, Array{Cell, 1}(), sym)
        
        cell
    end
end

function add_parent(cur::Cell, parent::Cell)
    if cur.parent == nothing
        cur.parent = parent
    else
        #insert given parent into chain
        parent.parent = cur.parent
        cur.parent = parent
    end
    push!(parent.children, cur)
end

function add_child(cur::Cell, child::Cell)
    push!(cur.children, child)
end

function show(io::IO, cell::Cell, ilevel::Int64=0)
    parent_present = cell.parent == nothing ? "none" : "present"
    sym_desc = cell.sym == nothing ? "(nothing)" : cell.sym

    iprintln(io, "Cell:", ilevel)

    iprintln(io, "energy: $(cell.energy)", ilevel + 1)
    iprintln(io, "parent: $(parent_present)", ilevel + 1)
    iprintln(io, "children: $(length(cell.children))", ilevel + 1)
    iprintln(io, "sym: $(sym_desc)", ilevel + 1)
    
    iprintln(io, "gene_states:", ilevel + 1)
    map(gs -> GeneStateMod.show(io, gs, ilevel + 1), cell.gene_states)
end

end

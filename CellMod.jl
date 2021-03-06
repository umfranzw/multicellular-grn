module CellMod

using RunMod
using GeneStateMod
using GeneMod
using ProteinStoreMod
using ProteinMod
using SymMod
using SymProbsMod
using MiscUtilsMod
using ProteinPropsMod
using Printf

import Base.show

export Cell

mutable struct Cell
    config::Config
    gene_states::Array{GeneState, 1}
    proteins::ProteinStore
    parent::Union{Cell, Nothing}
    children::Array{Cell, 1}
    probs::SymProbs
    sym::Union{Sym, Nothing}
    age::Int64
    id::Union{UInt64, Nothing} #for use in preventing self-binding, and in the UI

    function Cell(config::Config, genes::Array{Gene, 1}, age::Int64=0)
        gene_states = map(g -> GeneState(config.run, g), genes)
        proteins = ProteinStore(config.run)
        children = Array{Cell, 1}()

        cell = new(config, gene_states, proteins, nothing, children, SymProbs(), nothing, age, nothing)
        cell.id = hash(cell)
        
        cell
    end

    function Cell(config::Config, gene_states::Array{GeneState, 1}, age::Int64=0)
        cell = new(config, gene_states, ProteinStore(config.run), nothing, Array{Cell, 1}(), SymProbs(), nothing, age, nothing)
        
        cell
    end
end

function insert_initial_proteins(cell::Cell, proteins::Array{Protein, 1})
    i = 1
    while i <= min(length(proteins), cell.config.run.max_proteins_per_cell, cell.config.run.max_initial_proteins)
        #it is possible that not all initial proteins in the array are unique. That's ok, since they'll all (including the duplicates) be subject to evolution.
        #the ProteinStore's insert() method ensures that only one (the last one) protein in each pair of duplicates gets inserted.
        #note: this means some individual's root cells may have fewer initial proteins than others...
        #note: we push a copy so the indiv's initial_cell_proteins array stays intact as the simulation modifies protein's concs
        #in the root cell
        copy = deepcopy(proteins[i])
        copy.src_cell_id = cell.id
        ProteinStoreMod.insert(cell.proteins, copy)

        i += 1
    end
end

function add_parent(cur::Cell, parent::Cell)
    cur.parent = parent
    push!(parent.children, cur)
end

function add_child(cur::Cell, child::Cell)
    push!(cur.children, child)
end

function fix_sym(cell::Cell)
    cell.sym = SymProbsMod.choose_sym(cell.probs, cell.config)
end

function show(io::IO, cell::Cell, ilevel::Int64=0)
    parent_present = cell.parent == nothing ? "none" : "present"
    sym_desc = cell.sym == nothing ? "(nothing)" : cell.sym

    iprintln(io, "Cell:", ilevel)

    iprintln(io, "id: $(cell.id)", ilevel + 1)
    iprintln(io, "parent: $(parent_present)", ilevel + 1)
    iprintln(io, "num children: $(length(cell.children))", ilevel + 1)
    iprintln(io, "sym: $(sym_desc)", ilevel + 1)
    iprintln(io, "age: $(cell.age)", ilevel + 1)
    
    iprintln(io, "gene_states:", ilevel + 1)
    foreach(gs -> GeneStateMod.show(io, gs, ilevel + 1), cell.gene_states)

    SymProbsMod.show(io, cell.probs, ilevel + 2)
end

end

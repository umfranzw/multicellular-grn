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
    sensors::Dict{ProteinPropsMod.ProteinLoc, Array{Float64, 1}}
    sym::Union{Sym, Nothing}
    age::Int64
    id::Union{UInt64, Nothing}

    function Cell(config::Config, genes::Array{Gene, 1})
        gene_states = map(g -> GeneState(config, g), genes)
        proteins = ProteinStore()
        children = Array{Cell, 1}()

        cell = new(config, gene_states, proteins, nothing, children, SymProbs(), build_sensors(config.run, length(genes)), nothing, 0, nothing)
        cell.id = hash(cell)
        
        cell
    end

    function Cell(config::Config, gene_states::Array{GeneState, 1})
        cell = new(config, gene_states, ProteinStore(), nothing, Array{Cell, 1}(), SymProbs(), build_sensors(length(genes)), nothing, 0)
        
        cell
    end
end

function build_sensors(run::Run, num_concs::Int64)
    sensors = Dict{ProteinPropsMod.ProteinLoc, Array{Float64, 1}}()
    for loc in instances(ProteinPropsMod.ProteinLoc)
        sensors[loc] = repeat([run.initial_cell_sensor_conc], num_concs)
    end

    sensors
end

function adjust_sensor(cell::Cell, loc::ProteinPropsMod.ProteinLoc, delta::Array{Float64, 1})
    cell.sensors[loc] = clamp.(cell.sensors[loc] .+ delta, 0.0, cell.config.run.max_sensor_amount)
end

function insert_initial_proteins(cell::Cell, proteins::Array{Protein, 1})
    i = 1
    while i <= min(length(proteins), cell.config.run.max_proteins_per_cell)
        #it is possible that not all initial proteins in the array are unique. That's ok, since they'll all (including the duplicates) be subject to evolution.
        #the ProteinStore's insert() method ensures that only one (the last one) protein in each pair of duplicates gets inserted.
        #note: this means some individual's root cells may have fewer initial proteins than others...
        #note: we push a copy so the indiv's initial_cell_proteins array stays intact as the simulation modifies protein's concs
        #in the root cell
        copy = deepcopy(proteins[i])
        ProteinStoreMod.insert(cell.proteins, copy)

        i += 1
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

    iprintln(io, "sensors:", ilevel + 1)
    for loc in instances(ProteinPropsMod.ProteinLoc)
        name = string(loc)
        val = join(map(c -> @sprintf("%0.2f", c), cell.sensors[loc]), ", ")
        iprintln(io, "$(name): [$(val)]", ilevel + 2)
    end
    
    SymProbsMod.show(io, cell.probs, ilevel + 2)
end

end

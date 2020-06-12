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
    sensors::Dict{Int8, Array{Float64, 1}}
    sym::Union{Sym, Nothing}
    age::Int64
    id::Union{UInt64, Nothing} #for use in preventing self-binding, and in the UI

    function Cell(config::Config, genes::Array{Gene, 1}, age::Int64=0)
        gene_states = map(g -> GeneState(config.run, g), genes)
        proteins = ProteinStore(config.run)
        children = Array{Cell, 1}()

        cell = new(config, gene_states, proteins, nothing, children, SymProbs(), build_sensors(config.run, length(genes)), nothing, age, nothing)
        cell.id = hash(cell)
        
        cell
    end

    function Cell(config::Config, gene_states::Array{GeneState, 1}, age::Int64=0)
        cell = new(config, gene_states, ProteinStore(config.run), nothing, Array{Cell, 1}(), SymProbs(), build_sensors(length(genes)), nothing, age, nothing)
        
        cell
    end
end

function build_sensors(run::Run, num_concs::Int64)
    sensors = Dict{Int8, Array{Float64, 1}}()
    for loc in 0 : 2 + run.max_children #right, top, left, children
        sensors[Int8(loc)] = repeat([run.initial_cell_sensor_conc], num_concs)
    end

    sensors
end

function adjust_sensor(cell::Cell, loc::Int8, delta::Array{Float64, 1})
    #println(length(cell.sensors[loc]))
    #println(length(delta))
    cell.sensors[loc] = clamp.(cell.sensors[loc] .+ delta, 0.0, cell.config.run.max_sensor_amount)
end

function extend_sensors(cell::Cell, index::Int64)
    for loc in keys(cell.sensors)
        insert!(cell.sensors[loc], index, cell.config.run.initial_cell_sensor_conc)
    end
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

    iprintln(io, "sensors:", ilevel + 1)
    for loc in sort(collect(keys(cell.sensors)))
        name = string(loc)
        val = join(map(c -> @sprintf("%0.2f", c), cell.sensors[loc]), ", ")
        iprintln(io, "$(name): [$(val)]", ilevel + 2)
    end
    
    SymProbsMod.show(io, cell.probs, ilevel + 2)
end

end

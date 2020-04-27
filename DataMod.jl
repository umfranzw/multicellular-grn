module DataMod

using IndividualMod
using CellTreeMod
using CellMod
using RunMod
using CacheMod
using ProteinMod
using ProteinPropsMod
using ProteinStoreMod
using GeneMod
using SymMod
using ChainGraphMod
using GeneStateMod
using Statistics
using Printf
import TrackerMod
import Serialization
import CodecZlib

import Base.close

export Data

#tree_cache_size = 20
indiv_cache_size = 20

mutable struct Data
    #dictionary order is ea_step, index, reg_step
    #trees::Dict{Tuple{Int64, Int64, Int64}, CellTree}
    indivs::Cache{Tuple{Int64, Int64, Int64}, Individual}
    #(ea_step, index, reg_step) => (position, size)
    indivs_index::Dict{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64}}
    #dictionary order is ea_step, index
    #indivs::Dict{Tuple{Int64, Int64}, Individual}
    #indivs::Cache{Tuple{Int64, Int64}, Individual}
    #(ea_step, index) => (position, size)
    #indivs_index::Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}
    file_handle::IOStream
    run::Union{Run, Nothing}

    function Data(filename::String)
        file_path = join((RunMod.DATA_PATH, filename), "/")
        file_handle = open(file_path, "r")
        
        #ea_step, index, reg_step
        #trees = Dict{Tuple{Int64, Int64, Int64}, CellTree}()
        indivs = Cache{Tuple{Int64, Int64, Int64}, Individual}(DataMod.indiv_cache_size)
        indivs_index = Dict{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64}}()
        
        #ea_step, index
        #indivs = Dict{Tuple{Int64, Int64}, Individual}()
        #indivs = Cache{Tuple{Int64, Int64}, Individual}(DataMod.indiv_cache_size)
        #indivs_index = Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}()
        
        data = new(indivs, indivs_index, file_handle, nothing)
        data.run = get_run(data)
        create_index(data)

        data
    end
end

function close(data::Data)
    close(data.file_handle)
end

function get_run(data::Data)
    seek(data.file_handle, 0)
    tag_type = read(data.file_handle, TrackerMod.TagType)
    
    if tag_type != TrackerMod.RunState
        println("Error - the run is not the first thing in the data file.")
        close(data)
        exit(1)
    end
    
    size = read(data.file_handle, Int64)

    read_obj(data, position(data.file_handle), size)
end

# function get_indiv(data::Data, ea_step::Int64, pop_index::Int64)
#     key = (ea_step,  pop_index)
#     get_obj(data, key, :indivs, :indivs_index)
# end

function get_indiv(data::Data, ea_step::Int64, pop_index::Int64, reg_step::Int64)
    key = (ea_step, pop_index, reg_step)
    get_obj(data, key, :indivs, :indivs_index)
end

function get_obj(data::Data, key::Any, store_field::Symbol, index_field::Symbol)
    store = getfield(data, store_field)
    index = getfield(data, index_field)
    
    if key ∉ keys(store)
        pos, size = index[key]
        store[key] = read_obj(data, pos, size)
    end    
    
    store[key]
end

function read_obj(data::Data, pos::Int64, size::Int64)
    seek(data.file_handle, pos)

    comp_obj = read(data.file_handle, size)
    ser_obj = CodecZlib.transcode(CodecZlib.GzipDecompressor, comp_obj)
    ser_buf = IOBuffer(ser_obj)

    Serialization.deserialize(ser_buf)
end

function create_index(data::Data)
    while !eof(data.file_handle)
        tag_type = read(data.file_handle, TrackerMod.TagType)
        if tag_type == TrackerMod.RunState
            #note: no tag here
            size = read(data.file_handle, Int64)
            #no need to add this to an index - it's always at the start of the data file
            
        elseif tag_type == TrackerMod.IndivState
            ea_step = read(data.file_handle, Int64)
            reg_step = read(data.file_handle, Int64)
            pop_index = read(data.file_handle, Int64)
            size = read(data.file_handle, Int64)

            key = (ea_step, pop_index, reg_step)
            data.indivs_index[key] = (position(data.file_handle), size)
        end

        #go to next item
        seek(data.file_handle, position(data.file_handle) + size)
    end
end

#analysis functions

function get_conc_sum_for_tree(tree::CellTree, protein::Protein)
    get_concs_for_tree(tree, protein, sum)
end

function get_conc_max_for_tree(tree::CellTree, protein::Protein)
    get_concs_for_tree(tree, protein, maximum)
end

function get_conc_mean_for_tree(tree::CellTree, protein::Protein)
    get_concs_for_tree(tree, protein, mean)
end

function get_conc_excess_for_tree(tree::CellTree, protein::Protein, threshold::Float64)
    get_concs_for_tree(tree, protein, concs -> clamp(concs .- threshold, 0.0, 1.0))
end

function get_concs_for_tree(tree::CellTree, protein::Protein, fcn::Union{Function, Nothing}=nothing)
    result = Dict{Cell, Array{Float64, 1}}()
    CellTreeMod.traverse(cell -> result[cell] = get_concs_for_cell(cell, protein, fcn), tree)

    result
end

function get_concs_for_cell(cell::Cell, protein::Protein, fcn::Union{Function, Nothing}=nothing)
    result = nothing
    target = ProteinStoreMod.get(cell.proteins, protein.props)
    if target != nothing
        if fcn == nothing
            result = target.concs
        else
            result = [fcn(target.concs)]
        end
    end
    
    result    
end

function get_protein_info_for_indiv(indiv::Individual)
    info = Set{Tuple{String, String, String, String, UInt8, Bool, ProteinProps, UInt64}}()
    if indiv.cell_tree.root != nothing
        CellTreeMod.traverse(cell -> union!(info, get_protein_info_for_cell(cell)), indiv.cell_tree.root)
    end

    array = collect(info)
    sort!(array, by=row -> row[end])
    
    array
end

function get_protein_info_for_cell(cell::Cell)
    info = Set{Tuple{String, String, String, String, UInt8, Bool, ProteinProps, UInt64}}()
    proteins = ProteinStoreMod.get_all(cell.proteins)
    for protein in proteins
        item  = (
            string(protein.props.type),
            string(protein.props.fcn),
            string(protein.props.action),
            string(protein.props.loc),
            protein.props.arg,
            protein.is_initial,
            protein.props,
            hash((protein.props.type,
                  protein.props.fcn,
                  protein.props.action,
                  protein.props.loc,
                  protein.props.arg,
                  protein.is_initial)) #note: if the same protein is present in two different cells, they'll hash the same
            #this means they'll be considered identical in the gui
        )
        push!(info, item)
    end

    info
end

function get_probs_info_for_cell(cell::Cell)
    info = Array{Tuple{String, Float64}, 1}() #[(sym_label, prob), ...]
    for (sym, val) in cell.probs.probs
        label = SymMod.to_str(sym)
        push!(info, (label, val))
    end

    sort!(info, by=pair -> pair[1]) #sort by label
    
    info
end

function get_sensor_concs(cell::Cell)
    concs = Dict{String, Array{Float64, 1}}()
    for loc in instances(ProteinPropsMod.ProteinLoc)
        concs[string(loc)] = cell.sensors[loc]
    end

    concs
end

function build_graph_for_cell(data::Data, ea_step::Int64, pop_index::Int64, cell::Cell)
    graph = ChainGraph()

    for reg_step in 1:data.run.reg_steps + 1
        indiv = DataMod.get_indiv(data, ea_step, pop_index, reg_step)
        tree = indiv.cell_tree
        cur_cell = CellTreeMod.find_by_id(tree, cell.id)
        if cur_cell != nothing
            for gs in cur_cell.gene_states
                #find all bindings and insert them into the graph
                for bound_protein in gs.bindings
                    if bound_protein != nothing
                        ChainGraphMod.add_node(graph, (bound_protein.props, bound_protein.is_initial))
                        ChainGraphMod.add_node(graph, gs.gene)
                        ChainGraphMod.add_edge(graph, (bound_protein.props, bound_protein.is_initial), gs.gene)
                    end
                end

                #recalc all productions and insert them into the graph
                rates = GeneStateMod.get_prod_rates(gs)
                for (prod_index, rate) in rates
                    if rate > 0
                        prod_props = gs.gene.prod_sites[prod_index]
                        #the produced protein should have been inserted into the store already
                        #protein = ProteinStoreMod.get(cur_cell.proteins, prod_props)
                        #note: binding must occur in order for production to occur, so here, the gene is already in the graph
                        #(it was inserted above)
                        #Just need to add the protein and an edge
                        ChainGraphMod.add_node(graph, (prod_props, false)) #produced proteins are never initial
                        ChainGraphMod.add_edge(graph, gs.gene, (prod_props, false))
                    end
                end
            end
        end
    end

    ChainGraphMod.plot(graph)
end

function get_gs_table_data(data::Data, cell::Cell, ea_step::Int64, indiv_index::Int64, reg_step::Int64)
    #table data should be a 2d array
    #each row is a binding site
    #cols are as described in vis/CellArea.py/build_gene_states_tab
    table = Array{Array{String, 1}, 1}()

    for gene_index in 1:length(cell.gene_states)
        gs = cell.gene_states[gene_index]
        gene_index_str = string(gene_index)
        bind_logic = string(gs.gene.bind_logic)

        bind_sites = Array{String, 1}()
        thresholds = Array{String, 1}()
        bound_proteins = Array{String, 1}()
        for i in 1:length(gs.gene.bind_sites)
            push!(bind_sites, GeneMod.get_bind_site_str(gs.gene, i))
            push!(thresholds, @sprintf("%0.2f", gs.gene.bind_sites[i].threshold))
            if gs.bindings[i] == nothing
                push!(bound_proteins, "")
            else
                push!(bound_proteins, ProteinPropsMod.to_str(gs.bindings[i].props))
            end
        end

        prod_rates = GeneStateMod.get_prod_rates(cell.gene_states[gene_index])
        prod_sites = repeat([""], length(gs.gene.bind_sites))
        prod_rates_strs = repeat([""], length(gs.gene.bind_sites))
        for (prod_index, rate) in prod_rates
            prod_sites[prod_index] = GeneMod.get_prod_site_str(gs.gene, prod_index)
            prod_rates_strs[prod_index] = @sprintf("%0.2f", rate)
        end

        row = Array{String, 1}()
        push!(row, gene_index_str)
        push!(row, bind_logic)
        for i in 1:length(gs.gene.bind_sites)
            push!(row, bind_sites[i])
            push!(row, bound_proteins[i])
            push!(row, thresholds[i])
        end

        for i in 1:length(gs.gene.bind_sites)
            push!(row, prod_sites[i])
            push!(row, prod_rates_strs[i])
        end

        push!(table, row)
    end

    table
end

function get_best_fitnesses(data::Data)
    bests = Array{Float64, 1}()
    gen_bests = Array{Float64, 1}()
    gen_avgs = Array{Float64, 1}()
    final_reg_step = data.run.reg_steps + 1

    best = nothing
    for ea_step in 0 : data.run.step_range.step : data.run.ea_steps
        gen_avg = 0.0
        gen_best = nothing
        for pop_index in 1:data.run.pop_size
            indiv = DataMod.get_indiv(data, ea_step, pop_index, final_reg_step)
            gen_avg += indiv.fitness
            
            if gen_best == nothing || indiv.fitness < gen_best
                gen_best = indiv.fitness
                
                #note: only need to update overall best if gen_best is updated
                if best == nothing || gen_best < best
                    best = gen_best
                end
            end
        end
        gen_avg /= data.run.pop_size
        
        push!(bests, best)
        push!(gen_bests, gen_best)
        push!(gen_avgs, gen_avg)
    end

    (bests, gen_bests, gen_avgs)
end

end

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

cache_size = 10

mutable struct Data
    #dictionary order is ea_step, index, reg_step
    #trees::Dict{Tuple{Int64, Int64, Int64}, CellTree}
    trees::Cache{Tuple{Int64, Int64, Int64}, CellTree}
    #(ea_step, index, reg_step) => (position, size)
    trees_index::Dict{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64}}
    #dictionary order is ea_step, index
    #indivs::Dict{Tuple{Int64, Int64}, Individual}
    indivs::Cache{Tuple{Int64, Int64}, Individual}
    #(ea_step, index) => (position, size)
    indivs_index::Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}
    file_handle::IOStream
    run::Union{Run, Nothing}

    function Data(filename::String)
        file_path = join((RunMod.DATA_PATH, filename), "/")
        file_handle = open(file_path, "r")
        
        #ea_step, index, reg_step
        #trees = Dict{Tuple{Int64, Int64, Int64}, CellTree}()
        trees = Cache{Tuple{Int64, Int64, Int64}, CellTree}(DataMod.cache_size)
        trees_index = Dict{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64}}()
        
        #ea_step, index
        #indivs = Dict{Tuple{Int64, Int64}, Individual}()
        indivs = Cache{Tuple{Int64, Int64}, Individual}(DataMod.cache_size)
        indivs_index = Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}()
        
        data = new(trees, trees_index, indivs, indivs_index, file_handle, nothing)
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

function get_indiv(data::Data, ea_step::Int64, pop_index::Int64)
    key = (ea_step,  pop_index)
    get_obj(data, key, :indivs, :indivs_index)
end

function get_tree(data::Data, ea_step::Int64, pop_index::Int64, reg_step::Int64)
    key = (ea_step, pop_index, reg_step)
    get_obj(data, key, :trees, :trees_index)
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
            pop_index = read(data.file_handle, Int64)
            size = read(data.file_handle, Int64)
            
            key = (ea_step, pop_index)
            data.indivs_index[key] = (position(data.file_handle), size)

        elseif tag_type == TrackerMod.CellTreeState
            ea_step = read(data.file_handle, Int64)
            reg_step = read(data.file_handle, Int64)
            pop_index = read(data.file_handle, Int64)
            size = read(data.file_handle, Int64)

            key = (ea_step, pop_index, reg_step)
            data.trees_index[key] = (position(data.file_handle), size)
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

function get_protein_info_for_tree(tree::CellTree)
    info = Set{Tuple{String, String, String, String, Int8, Bool, ProteinProps, UInt64}}()
    if tree.root != nothing
        CellTreeMod.traverse(cell -> union!(info, get_protein_info_for_cell(cell)), tree.root)
    end

    array = collect(info)
    sort!(array, by=row -> row[end])
    
    array
end

function get_protein_info_for_cell(cell::Cell)
    info = Set{Tuple{String, String, String, String, Int8, Bool, ProteinProps, UInt64}}()
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
        tree = DataMod.get_tree(data, ea_step, pop_index, reg_step)
        cur = CellTreeMod.find_by_id(tree, cell.id)
        if cur != nothing
            for gs in cur.gene_states
                #find all bindings and insert them into the graph
                for bound_protein in gs.bindings
                    if bound_protein != nothing
                        ChainGraphMod.add_node(graph, bound_protein)
                        ChainGraphMod.add_node(graph, gs.gene)
                        ChainGraphMod.add_edge(graph, bound_protein, gs.gene)
                    end
                end

                #find all productions and insert them into the graph
                rates = GeneStateMod.get_prod_rates(gs)
                for (prod_index, rate) in rates
                    if rate > 0
                        prod_props = gene.prod_sites[prod_index]
                        #the produced protein should have been inserted into the store already
                        protein = ProteinStoreMod.get(cur.proteins, prod_props)
                        #note: binding must occur in order for production to occur, so here, the gene is already in the graph
                        #(it was inserted above)
                        #Just need to add the protein and an edge
                        ChainGraphMod.add_node(graph, protein)
                        ChainGraphMod.add_edge(gs.gene, protein)
                    end
                end
            end
        end
    end

    #println(ChainGraphMod.gen_dot_code(graph))
    ChainGraphMod.plot(graph)
end

function get_cell_from_prev_reg_step(data::Data, cur_cell::Cell, ea_step::Int64, indiv_index::Int64, cur_reg_step::Int64)
    prev_cell = nothing

    if cur_reg_step != 0
        prev_tree = get_tree(data, ea_step, indiv_index, cur_reg_step - 1)
        prev_cell = CellTreeMod.find_by_id(prev_tree, cur_cell.id)
    end

    prev_cell #will be nothing if not present
end

function get_gs_table_data(data::Data, cell::Cell, ea_step::Int64, indiv_index::Int64, reg_step::Int64)
    #table data should be a 2d array
    #each row is a binding site
    #cols are as described in vis/CellArea.py/build_gene_states_tab
    table = Array{Array{String, 1}, 1}()

    prev_cell = DataMod.get_cell_from_prev_reg_step(data, cell, ea_step, indiv_index, reg_step)
    for gene_index in 1:length(cell.gene_states)
        gs = cell.gene_states[gene_index]
        gene_index_str = string(gene_index)
        bind_logic = string(gs.gene.bind_logic)

        bind_sites = Array{String, 1}()
        bound_proteins = Array{String, 1}()
        for i in 1:length(gs.gene.bind_sites)
            push!(bind_sites, GeneMod.get_bind_site_str(gs.gene, i))
            if gs.bindings[i] == nothing
                push!(bound_proteins, "")
            else
                push!(bound_proteins, ProteinPropsMod.to_str(gs.bindings[i].props))
            end
        end

        if prev_cell == nothing
            prod_rates = Array{NamedTuple{(:prod_index, :rate), Tuple{Int64, Float64}}, 1}()
        else
            prod_rates = GeneStateMod.get_prod_rates(prev_cell.gene_states[gene_index])
        end
        prod_sites = repeat([""], length(gs.gene.bind_sites))
        prod_rates_strs = repeat([""], length(gs.gene.bind_sites))
        for (prod_index, rate) in prod_rates
            prod_sites[prod_index] = GeneMod.get_prod_site_str(gs.gene, prod_index)
            prod_rates_strs[prod_index] = @printf("%0.2f", rate)
        end

        row = Array{String, 1}()
        push!(row, gene_index_str)
        push!(row, bind_logic)
        for i in 1:length(gs.gene.bind_sites)
            push!(row, bind_sites[i])
            push!(row, bound_proteins[i])
        end

        for i in 1:length(gs.gene.bind_sites)
            push!(row, prod_sites[i])
            push!(row, prod_rates_strs[i])
        end

        push!(table, row)
    end

    table
end

end

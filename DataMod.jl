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
using InternalGraphMod
using NeighbourGraphMod
using GeneStateMod
using Statistics
using Printf
using TrackerMod
using Formatting
import Serialization
import CodecZlib

import Base.close

export Data

indiv_cache_size = 20

mutable struct Data
    #dictionary order is ea_step, index, reg_step
    indivs::Cache{Tuple{Int64, Int64, Int64}, Individual}
    #(ea_step, index, reg_step) => (position, size)
    indivs_index::Dict{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64}}
    #ea_step, pop_index
    fitnesses::Union{Array{Array{Float64, 1}, 1}, Nothing}
    file_handle::IOStream
    run_best::Union{BestInfo, Nothing}
    run::Union{Run, Nothing}

    function Data(filename::String)
        file_path = join((RunMod.DATA_PATH, filename), "/")
        file_handle = open(file_path, "r")
        
        indivs = Cache{Tuple{Int64, Int64, Int64}, Individual}(DataMod.indiv_cache_size)
        indivs_index = Dict{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64}}()
        
        
        data = new(indivs, indivs_index, nothing, file_handle, nothing, nothing)
        #data.run = get_run(data)
        create_index(data) #note: this also initializes run_best and run

        data
    end
end

function close(data::Data)
    close(data.file_handle)
end

# function get_run(data::Data)
#     seek(data.file_handle, 0)
#     tag_type = read(data.file_handle, TrackerMod.TagType)

#     if tag_type != TrackerMod.RunState
#         println("Error - the run is not the first thing in the data file.")
#         close(data)
#         exit(1)
#     end

#     size = read(data.file_handle, Int64)

#     read_obj(data, position(data.file_handle), size)
# end

# function get_indiv(data::Data, ea_step::Int64, pop_index::Int64)
#     key = (ea_step,  pop_index)
#     get_obj(data, key, :indivs, :indivs_index)
# end

function get_indiv(data::Data, ea_step::Int64, pop_index::Int64, reg_step::Int64)
    key = (ea_step, pop_index, reg_step)
    if key == data.run_best.index
        indiv = data.run_best.indiv
    else
        indiv = get_obj(data, key)
    end

    indiv
end

function get_obj(data::Data, key::Any)
    if key ∉ keys(data.indivs)
        pos, size = data.indivs_index[key]
        data.indivs[key] = read_obj(data, pos, size)
    end    
    
    data.indivs[key]
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
            #no need to add this to an index - we'll save it in the data object
            data.run = read_obj(data, position(data.file_handle), size)
            
        elseif tag_type == TrackerMod.IndivState
            ea_step = read(data.file_handle, Int64)
            reg_step = read(data.file_handle, Int64)
            pop_index = read(data.file_handle, Int64)
            size = read(data.file_handle, Int64)

            key = (ea_step, pop_index, reg_step)
            data.indivs_index[key] = (position(data.file_handle), size)

            #go to next item
            seek(data.file_handle, position(data.file_handle) + size)

        elseif tag_type == TrackerMod.RunBestInfoState
            #note: no tag here
            size = read(data.file_handle, Int64)
            #no need to add this to an index - we'll save it in the data object
            data.run_best = read_obj(data, position(data.file_handle), size)

        elseif tag_type == TrackerMod.FitnessesState
            #note: no tag here
            size = read(data.file_handle, Int64)
            #no need to add this to an index - we'll save it in the data object
            data.fitnesses = read_obj(data, position(data.file_handle), size)
        end
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
    info = Dict{ProteinProps, Tuple{String, String, String, Int8, Bool, ProteinProps, UInt64}}()
    if indiv.cell_tree.root != nothing
        CellTreeMod.traverse(cell -> merge!(info, get_protein_info_for_cell(cell)), indiv.cell_tree.root)
    end

    array = collect(values(info))
    sort!(array, by=row -> row[end])
    
    array
end

function get_protein_info_for_cell(cell::Cell)
    info = Dict{ProteinProps, Tuple{String, String, String, Int8, Bool, ProteinProps, UInt64}}()
    proteins = ProteinStoreMod.get_all(cell.proteins)
    for protein in proteins
        item  = (
            string(protein.props.type),
            string(protein.props.tag),
            string(protein.props.action),
            protein.props.arg,
            protein.is_initial,
            protein.props,
            hash(protein.props) #note: if the same protein is present in two different cells, they'll hash the same
            #this means they'll be considered identical in the gui
        )
        info[protein.props] = item
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
    concs = Array{Array{Float64, 1}, 1}()
    num_locs = 3 + cell.config.run.max_children
    for loc in 0 : num_locs - 1
        push!(concs, cell.sensors[loc])
    end

    concs
end

function add_neighbour_edges(graph::NeighbourGraph, cell::Cell, id_to_cell::Dict{UInt64, Cell})
    for gs in cell.gene_states
        for bound_protein in gs.bindings
            if bound_protein != nothing && bound_protein.props.type == ProteinPropsMod.Neighbour
                #note: this assumes there are no initial neighbour proteins (there's always a src cell)
                src = id_to_cell[bound_protein.src_cell_id]
                label = ProteinPropsMod.to_str(bound_protein.props, bound_protein.is_initial)
                NeighbourGraphMod.add_edge(graph, src, cell)
            end
        end
    end
end

function build_neighbour_comm_graph(data::Data, ea_step::Int64, pop_index::Int64, reg_step::Int64)
    indiv = DataMod.get_indiv(data, ea_step, pop_index, reg_step)
    graph = NeighbourGraph(indiv.cell_tree)
    id_to_cell = Dict{UInt64, Cell}()
    CellTreeMod.traverse(cell -> id_to_cell[cell.id] = cell, indiv.cell_tree)
    CellTreeMod.traverse(cell -> add_neighbour_edges(graph, cell, id_to_cell), indiv.cell_tree)
    
    NeighbourGraphMod.plot(graph)
end

function build_graph_for_cell(data::Data, ea_step::Int64, pop_index::Int64, reg_step::Int64, cell::Cell)
    graph = InternalGraph()
    indiv = DataMod.get_indiv(data, ea_step, pop_index, reg_step)
    tree = indiv.cell_tree
    cur_cell = CellTreeMod.find_by_id(tree, cell.id)
    if cur_cell != nothing
        #add all genes to the graph
        for gs in cur_cell.gene_states
            InternalGraphMod.add_gene(graph, gs.gene)
        end
        
        for gs in cur_cell.gene_states
            #find all bindings and insert them into the graph
            for bound_protein in gs.bindings
                if bound_protein != nothing
                    InternalGraphMod.add_props(graph, bound_protein.props, bound_protein.is_initial)
                    InternalGraphMod.add_edge(graph, bound_protein.props, gs.gene)
                end
            end

            #recalc all productions and insert them into the graph
            rates = GeneStateMod.get_prod_rates(gs)
            for (prod_index, rate) in rates
                if rate > 0
                    prod_site = gs.gene.prod_sites[prod_index]
                    prod_props = ProteinProps(prod_site.type, prod_site.tag, prod_site.action, prod_site.arg)
                    #Just need to add the protein and an edge
                    InternalGraphMod.add_props(graph, prod_props, false) #note: produced proteins are never initial
                    InternalGraphMod.add_edge(graph, gs.gene, prod_props)
                end
            end
        end
    end

    InternalGraphMod.plot(graph)
end

function save_all_graphs_for_cell(data::Data, ea_step::Int64, pop_index::Int64, cell::Cell, path::String)
    if !endswith(path, '/')
        path = string(path, '/')
    end

    digits_needed = Int64(floor(log10(data.run.reg_steps + 1)))
    fmt_spec = Formatting.FormatSpec("0$(digits_needed)d")
    for reg_step in 1:data.run.reg_steps + 1
        png_data = build_graph_for_cell(data, ea_step, pop_index, reg_step, cell)
        if png_data != nothing
            filename = string(path, Formatting.fmt(fmt_spec, reg_step), ".png")
            handle = open(filename, "w")
            write(handle, png_data)
            close(handle)
        end
    end
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
                push!(bound_proteins, ProteinPropsMod.to_str(gs.bindings[i].props, gs.bindings[i].is_initial))
            end
        end

        prod_rates = GeneStateMod.get_prod_rates(cell.gene_states[gene_index])
        prod_sites = repeat([""], length(gs.gene.bind_sites))
        prod_rates_strs = repeat([""], length(gs.gene.bind_sites))
        for i in 1:length(gs.gene.bind_sites)
            prod_sites[i] = GeneMod.get_prod_site_str(gs.gene, i)
        end
        
        for (prod_index, rate) in prod_rates
            #prod_sites[prod_index] = GeneMod.get_prod_site_str(gs.gene, prod_index)
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

#saves for all reg steps
function save_all_gs_table_data(data::Data, cell::Cell, ea_step::Int64, indiv_index::Int64, filename::String)
    handle = open(filename, "w")
    headers = ["Reg Step", "Gene Index", "Bind Logic"]
    for i in 1:data.run.bind_sites_per_gene
        push!(headers, "Bind Site $(i)")
        push!(headers, "Bound Protein $(i)")
    end
    push!(headers, "Threshold")
    for i in 1:data.run.bind_sites_per_gene
        push!(headers, "Prod Site $(i)")
        push!(headers, "Prod Rate $(i)")
    end
    write(handle, join(headers, ","))
    write(handle, "\n")
    
    for reg_step in 1:data.run.reg_steps + 1
        indiv = DataMod.get_indiv(data, ea_step, indiv_index, reg_step)
        tree = indiv.cell_tree
        cur_cell = CellTreeMod.find_by_id(tree, cell.id)

        if cur_cell != nothing
            step_data = get_gs_table_data(data, cur_cell, ea_step, indiv_index, reg_step)
            for row in step_data
                write(handle, "$(reg_step),")
                write(handle, join(row, ","))
                write(handle, "\n")
            end
        end
    end

    close(handle)
end

function get_best_fitnesses(data::Data)
    bests = Array{Float64, 1}()
    gen_bests = Array{Float64, 1}()
    gen_avgs = Array{Float64, 1}()

    best = nothing
    for ea_step in 1:length(data.fitnesses)
        gen_avg = 0.0
        gen_best = nothing
        for pop_index in 1:data.run.pop_size
            fitness = data.fitnesses[ea_step][pop_index]
            gen_avg += fitness

            if gen_best == nothing || fitness < gen_best
                gen_best = fitness

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

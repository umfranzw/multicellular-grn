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
using BestInfoMod
using FitnessInfoMod
import CheckpointMod
import Serialization
import CompressionMod

import Base.close

export Data

indiv_cache_size = 25

struct IndexKey
    #ea_step, pop_index, reg_step
    step_tag::Tuple{Int64, Int64, Int64}
    state_time::Union{TrackerMod.StateTime, Nothing}
end

struct IndexInfo
    state_content_type::Union{TrackerMod.StateContentType, Nothing}
    size::Int64
    file_pos::Int64
end

mutable struct Data
    indivs::Cache{IndexKey, Individual}
    indivs_index::Dict{IndexKey, IndexInfo}
    fitnesses::Union{Array{Array{Float64, 1}, 1}, Nothing} #indexing: ea_step, pop_index
    file_handle::IOStream
    run_best::Union{BestInfo, Nothing}
    run::Union{Run, Nothing}

    function Data(filename::String)
        file_path = join((RunMod.DATA_PATH, filename), "/")
        file_handle = open(file_path, "r")

        #determine the compression method that was used on the file
        #this is the first thing in the file - it's of type UInt8 (one of the values from the CompressionMod.CompressionAlg enum)
        compression_alg = CompressionMod.CompressionAlg(read(file_handle, UInt8))
        
        indivs = Cache{IndexKey, Individual}(DataMod.indiv_cache_size)
        indivs_index = Dict{IndexKey, IndexInfo}()
        
        data = new(indivs, indivs_index, nothing, file_handle, nothing, nothing)
        create_index(data, compression_alg) #note: this also initializes run_best and run

        data
    end
end

function close(data::Data)
    close(data.file_handle)
end

function get_indiv(data::Data, ea_step::Int64, pop_index::Int64, reg_step::Int64, state_time::TrackerMod.StateTime)
    get_indiv(data, IndexKey((ea_step, pop_index, reg_step), state_time))
end

function get_indiv(data::Data, key::IndexKey)
    if data.run_best != nothing && key.step_tag == data.run_best.index
        indiv = data.run_best.indiv
    else
        #if indiv is not in the cache, put it into the cache
        if key ∉ keys(data.indivs)
            index_info = data.indivs_index[key]
            #check if this is a checkpoint indiv
            #if so, the data stored in the file is the individual. Retrieve it.
            if index_info.state_content_type == TrackerMod.Checkpoint
                data.indivs[key] = read_obj(data, index_info.file_pos, index_info.size, data.run.compression_alg)
                #otherwise, we need to rebuild it using the last checkpoint
            else
                #the data stored in the file is a compressed CheckpointMod.ChangInfo struct. Decompress and deserialize it.
                change_info = read_obj(data, index_info.file_pos, index_info.size, data.run.compression_alg)
                #println("change_info: $(change_info)")
                
                #find the reg_step that the last checkpoint was recorded on
                ea_step, pop_index, reg_step = key.step_tag
                it_pair = iterate(TrackerMod.checkpoint_reg_steps)
                checkpt_reg_step = -1
                while checkpt_reg_step < reg_step && it_pair != nothing
                    checkpt_reg_step, state = it_pair
                    it_pair = iterate(TrackerMod.checkpoint_reg_steps, state)
                end

                if checkpt_reg_step == -1
                    throw(Error("The first checkpoint regstep must be 1"))
                end
                
                #get the index_info for the checkpoint indiv from data.indivs_index
                #note: only AfterBind indivs can be checkpoints - all AfterProd indivs use the last Afterbind checkpoint
                checkpt_key = IndexKey((ea_step, pop_index, checkpt_reg_step), TrackerMod.AfterBind)
                #println("checkpt_key: $(ea_step), $(pop_index), $(checkpt_reg_step)")
                checkpt_info = data.indivs_index[checkpt_key]
                #read the (compressed) bytes for the checkpoint indiv from the file
                #don't decompress here, since Checkpoint compression was performed on *compressed* data
                checkpt_bytes = read_bytes(data, checkpt_info.file_pos, checkpt_info.size, nothing)
                #reconstruct the target individual using the checkpoint bytes and the change_info
                resurrected_bytes = CheckpointMod.decompress(checkpt_bytes, change_info)
                decomp_resurrected_bytes = CompressionMod.decompress(data.run.compression_alg, resurrected_bytes)
                resurrected_indiv = Serialization.deserialize(IOBuffer(decomp_resurrected_bytes))

                #add it to the index
                data.indivs[key] = resurrected_indiv
            end
            
        end

        #retrieve the indiv from the cache
        indiv = data.indivs[key]
    end

    indiv
end

#will only perform decompression if decompression_alg is supplied
function read_obj(data::Data, pos::Int64, size::Int64, decompression_alg::Union{CompressionMod.CompressionAlg, Nothing}=nothing)
    decomp_bytes = read_bytes(data, pos, size, decompression_alg)
    buf = IOBuffer(decomp_bytes)

    Serialization.deserialize(buf)
end

#will only perform decompression if decompression_alg is supplied
function read_bytes(data::Data, pos::Int64, size::Int64, decompression_alg::Union{CompressionMod.CompressionAlg, Nothing}=nothing)
    seek(data.file_handle, pos)
    comp_bytes = read(data.file_handle, size)

    if decompression_alg == nothing
        result = comp_bytes
    else
        result = CompressionMod.decompress(decompression_alg, comp_bytes)
    end

    result
end

function create_index(data::Data, compression_alg::CompressionMod.CompressionAlg)
    #format: state_type, size, compressed_data, [step_tag], [state_time], [state_content_type]
    while !eof(data.file_handle)
        #read the information that is always present,
        #with the exception of the compressed data - instead, just record its position in the file
        state_type = read(data.file_handle, TrackerMod.StateType)
        size = read(data.file_handle, Int64)
        comp_data_pos = position(data.file_handle)
        
        if state_type == TrackerMod.RunState
            #no need to add this to an index - we'll save it in the data object
            data.run = read_obj(data, comp_data_pos, size, compression_alg)
            #println("read run")
            
        elseif state_type == TrackerMod.IndivState
            #since we don't want to consume the comp_data (just index it for later use), we need to skip over it to go to the following metadata
            seek(data.file_handle, position(data.file_handle) + size)
            
            #read step tag (ea_step, pop_index, reg_step)
            ea_step = read(data.file_handle, Int64)
            pop_index = read(data.file_handle, Int64)
            reg_step = read(data.file_handle, Int64)

            #read state_time and state_content_type
            state_time = TrackerMod.StateTime(read(data.file_handle, UInt8))
            state_content_type = read(data.file_handle, TrackerMod.StateContentType)

            key = IndexKey((ea_step, pop_index, reg_step), state_time)
            val = IndexInfo(state_content_type, size, comp_data_pos)
            data.indivs_index[key] = val
            #println("read indiv: ($ea_step, $pop_index, $reg_step)")

        elseif state_type == TrackerMod.RunBestInfoState
            #no need to add this to an index - we'll save it in the data object
            data.run_best = read_obj(data, comp_data_pos, size, compression_alg)
            #println("read run_best")

        elseif state_type == TrackerMod.FitnessesState
            #no need to add this to an index - we'll save it in the data object
            data.fitnesses = read_obj(data, comp_data_pos, size, compression_alg)
            #println("read fitnesses")
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
        for type in keys(gs.bindings)
            for bound_protein in gs.bindings[type]
                if bound_protein != nothing && bound_protein.props.type == ProteinPropsMod.Neighbour
                    #note: this assumes there are no initial neighbour proteins (there's always a src cell)
                    src = id_to_cell[bound_protein.src_cell_id]
                    label = ProteinPropsMod.to_str(bound_protein.props, bound_protein.is_initial)
                    NeighbourGraphMod.add_edge(graph, src, cell)
                end
            end
        end
    end
end

function build_neighbour_comm_graph(data::Data, ea_step::Int64, pop_index::Int64, reg_step::Int64)
    key = IndexKey((ea_step, pop_index, reg_step), TrackerMod.AfterBind)
    indiv = DataMod.get_indiv(data, key)
    graph = NeighbourGraph(indiv.cell_tree)
    id_to_cell = Dict{UInt64, Cell}()
    CellTreeMod.traverse(cell -> id_to_cell[cell.id] = cell, indiv.cell_tree)
    CellTreeMod.traverse(cell -> add_neighbour_edges(graph, cell, id_to_cell), indiv.cell_tree)
    
    NeighbourGraphMod.plot(graph)
end

function build_graph_for_cell(data::Data, ea_step::Int64, pop_index::Int64, reg_step::Int64, cell::Cell)
    graph = InternalGraph()
    key = IndexKey((ea_step, pop_index, reg_step), TrackerMod.AfterBind)
    indiv = DataMod.get_indiv(data, key)
    tree = indiv.cell_tree
    cur_cell = CellTreeMod.find_by_id(tree, cell.id)
    if cur_cell != nothing
        #add all genes to the graph
        for gs in cur_cell.gene_states
            InternalGraphMod.add_gene(graph, gs.gene)
        end
        
        for gs in cur_cell.gene_states
            #find all bindings and insert them into the graph
            for type in keys(gs.bindings)
                for bound_protein in gs.bindings[type]
                    if bound_protein != nothing
                        InternalGraphMod.add_props(graph, bound_protein.props, bound_protein.is_initial)
                        InternalGraphMod.add_edge(graph, bound_protein.props, gs.gene)
                    end
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

function export_table(data::Data, filename::String, data_fcn::Function)
    rows, max_cols = data_fcn()

    handle = open(filename, "w")
    for row in rows
        row_cols = length(row) - 1
        if row_cols < max_cols
            append!(row, repeat([","], max_cols - row_cols))
        end
        write(handle, join(row, ","))
        write(handle, "\n")
    end

    close(handle)
end

function export_gene_descs(data::Data, indiv_index::Int64, filename::String)
    export_table(data, filename, () -> get_gene_descs(data, indiv_index))
end

function export_reg_sim_info(data::Data, indiv_index::Int64, filename::String)
    export_table(data, filename, () -> get_reg_sim_info(data, indiv_index))
end

function export_fitness_info(data::Data, invid_index::Int64, filename::String)
    export_table(
        data,
        filename,
        function()
        table = get_fitness_info(data, indiv_index)
        cols = length(table[1]) #length of header row

        (table, cols)
        end
    )
end

function get_fitness_info(data::Data, pop_index::Int64)
    rows = Array{Array{String, 1}, 1}()
    headers = Array{String, 1}([
        "EA Step",
        "ID",
        "Fitness"
    ])
    for name in fieldnames(FitnessInfo)
        push!(headers, string(name))
    end
    push!(rows, headers)

    for ea_step in 0:data.run.ea_steps
        key = IndexKey((ea_step, pop_index, data.run.reg_steps + 1), TrackerMod.AfterBind)
        cur_indiv = DataMod.get_indiv(data, key)
        row = Array{String, 1}()
        #ea_step
        push!(row, string(ea_step))

        #ID
        push!(row, IndividualMod.get_id_str(cur_indiv))
        
        #fitness
        #push!(row, @sprintf("%0.2f", data.fitnesses[ea_step + 1][pop_index])) #note: since we store fitnesses for ea_step 0, the array is offset by one
        push!(row, @sprintf("%0.2f", cur_indiv.fitness))

        #fitness_info
        for name in fieldnames(FitnessInfo)
            push!(row, @sprintf("%0.2f", getproperty(cur_indiv.fitness_info, name)))
        end

        push!(rows, row)
    end
    
    rows
end

function get_reg_sim_info(data::Data, pop_index::Int64)
    #get reg_sim_info after the simulation has ended
    max_genes = -1
    rows = Array{Array{String, 1}, 1}()
    
    for ea_step in 0:data.run.ea_steps
        key = IndexKey((ea_step, pop_index, data.run.reg_steps + 1), TrackerMod.AfterBind)
        cur_indiv = DataMod.get_indiv(data, key)
        max_genes = max(max_genes, length(cur_indiv.genes))
        row = Array{String, 1}()
        #ea_step
        push!(row, string(ea_step))

        #ID
        push!(row, IndividualMod.get_id_str(cur_indiv))
        
        #fitness
        #push!(row, @sprintf("%0.2f", data.fitnesses[ea_step + 1][pop_index])) #note: since we store fitnesses for ea_step 0, the array is offset by one
        push!(row, @sprintf("%0.2f", cur_indiv.fitness))

        push!(row, string(cur_indiv.reg_sim_info.division_count))
        push!(row, string(cur_indiv.reg_sim_info.alter_sym_prob_count))
        
        #reg_sim_info
        for gene_index in 1:length(cur_indiv.genes)
            push!(row, string(cur_indiv.reg_sim_info.bind_count[gene_index]))
            push!(row, string(cur_indiv.reg_sim_info.produce_count[gene_index]))
        end

        push!(rows, row)
    end

    headers = Array{String, 1}()
    push!(headers, "EA Step")
    push!(headers, "ID")
    push!(headers, "Fitness")
    
    push!(headers, "Division Count")
    push!(headers, "Alter Sym Prob Count")
    
    for i in 1:max_genes
        push!(headers, "Bind Count $(i)")
        push!(headers, "Prod Count $(i)")
    end
    insert!(rows, 1, headers)

    (rows, max_genes)
end

function get_gene_descs(data::Data, pop_index::Int64)
    max_genes = -1
    rows = Array{Array{String, 1}, 1}()
    
    for ea_step in 0:data.run.ea_steps
        key = IndexKey((ea_step, pop_index, data.run.reg_steps + 1), TrackerMod.AfterBind)
        cur_indiv = DataMod.get_indiv(data, key)
        max_genes = max(max_genes, length(cur_indiv.genes))
        row = Array{String, 1}()
        #ea_step
        push!(row, string(ea_step))

        #ID
        push!(row, IndividualMod.get_id_str(cur_indiv))
        
        #fitness
        #push!(row, @sprintf("%0.2f", data.fitnesses[ea_step + 1][pop_index])) #note: since we store fitnesses for ea_step 0, the array is offset by one
        push!(row, @sprintf("%0.2f", cur_indiv.fitness))

        #gene descriptions
        for gene_index in 1:length(cur_indiv.genes)
            site_str = GeneMod.get_sites_str(cur_indiv.genes[gene_index])
            push!(row, site_str)
        end

        push!(rows, row)
    end

    headers = Array{String, 1}()
    push!(headers, "EA Step")
    push!(headers, "ID")
    push!(headers, "Fitness")
    
    for i in 1:max_genes
        push!(headers, "Gene $(i)")
    end
    insert!(rows, 1, headers)

    (rows, max_genes)
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
            if gs.bindings[BindSite][i] == nothing
                push!(bound_proteins, "")
            else
                push!(bound_proteins, ProteinPropsMod.to_str(gs.bindings[BindSite][i].props, gs.bindings[BindSite][i].is_initial))
            end
        end

        prod_rates = GeneStateMod.get_prod_rates(cell.gene_states[gene_index])
        prod_sites = repeat([""], length(gs.gene.bind_sites))
        prod_rates_strs = repeat([""], length(gs.gene.bind_sites))
        bound_inhib_proteins = Array{String, 1}()
        for i in 1:length(gs.gene.bind_sites)
            prod_sites[i] = GeneMod.get_prod_site_str(gs.gene, i)
            if gs.bindings[ProdSite][i] == nothing
                push!(bound_inhib_proteins, "")
            else
                push!(bound_inhib_proteins, ProteinPropsMod.to_str(gs.bindings[ProdSite][i].props, gs.bindings[ProdSite][i].is_initial))
            end
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
            push!(row, bound_inhib_proteins[i])
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
        key = IndexKey((ea_step, indiv_index, reg_step), TrackerMod.AfterBind)
        indiv = DataMod.get_indiv(data, key)
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

function get_indiv_fitness(data::Data, ea_step::Int64, pop_index::Int64)
    data.fitnesses[ea_step + 1][pop_index]
end

function get_best_fitnesses(data::Data)
    bests = Array{Float64, 1}()
    gen_bests = Array{Float64, 1}()
    gen_avgs = Array{Float64, 1}()

    best = nothing
    for ea_step in 0:data.run.ea_steps
        gen_avg = 0.0
        gen_best = nothing
        for pop_index in 1:data.run.pop_size
            fitness = data.fitnesses[ea_step + 1][pop_index] #remember, we store fitness for ea_step 0, so everything's offset
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

function get_base_seed(data::Data)
    key = IndexKey((0, 1, 1), TrackerMod.AfterBind)
    first_indiv = DataMod.get_indiv(data, key)
    
    first_indiv.config.rng.seed[1] - 1 #note: offset starts at 1
end

function get_pop_fitness_breakdown(data::Data)
    rows = Array{Array{String, 1}, 1}()
    sums = Dict{Symbol, Float64}()
    headers = Array{String, 1}()
    push!(headers, "EA Step")
    push!(headers, "Avg. Fitness")
    fitness_sum = 0
    for name in fieldnames(FitnessInfo)
        push!(headers, string(name))
    end
    
    push!(rows, headers)
    
    for ea_step in 0:data.run.ea_steps
        row = Array{String, 1}()
        fitness_sum = 0
        push!(row, string(ea_step))
        
        for pop_index in 1:data.run.pop_size
            cur_indiv = DataMod.get_indiv(data, ea_step, pop_index, data.run.reg_steps + 1, TrackerMod.AfterBind)
            fitness_sum += cur_indiv.fitness
            for name in fieldnames(FitnessInfo)
                if pop_index == 1
                    sums[name] = getproperty(cur_indiv.fitness_info, name)
                else
                    sums[name] += getproperty(cur_indiv.fitness_info, name)
                end
            end
        end

        push!(row, @sprintf("%0.2f", fitness_sum / data.run.pop_size))
        
        for name in fieldnames(FitnessInfo)
            push!(row, @sprintf("%0.2f", sums[name] / data.run.pop_size))
        end

        push!(rows, row)
    end

    rows
end

end

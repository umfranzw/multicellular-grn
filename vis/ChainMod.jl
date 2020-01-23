module ChainMod

using DataMod
using CellTreeMod
using ProteinPropsMod
using GeneMod
using CellMod

function find_chains(data::Data, ea_step::Int64, index::Int64, cell_index::Int64)
    links = find_links(data, ea_step, index, cell_index)
    chains = Array{Array{Union{ProteinProps, Gene}}, 1}()

    
end

function find_links(data::Data, ea_step::Int64, index::Int64, cell_index::Int64)
    DataMod.read_chunks_for_step(data, ea_step)

    #reg_step, gene_index, (gene, bindings, products)
    links = Array{Array{Tuple{Gene, Array{ProteinProps, 1}, Array{ProteinProps, 1}}}, 1}()

    for reg_step in 1:data.run.reg_steps
        tree = data.trees[ea_step][index][reg_step]
        cell = CellTreeMod.get_bf_node(tree, cell_index)
        if cell != nothing
            step_links = find_chains(cell)
            push!(links, step_links)
        end
    end

    links
end

function find_links(cell::Cell, graph::ChainGraph)
    #gene_index, (gene, bindings, products)
    links = Array{Tuple{Gene, Array{ProteinProps, 1}, Array{ProteinProps, 1}}}()
    for gene_index in 1:length(cell.gene_states)
        gs = cell.gene_states[gene_index]
        bindings = Array{ProteinProps, 1}()
        for protein in gs.reg_site_bindings
            if protein != nothing
                push!(bindings, protein.props)
            end
        end

        products = Array{ProteinProps, 1}()
        prod_rates = GeneStateMod.get_prod_rates(gs)
        for (rate, props) in zip(prod_rates, gs.prod_sites)
            if rate != nothing
                push!(products, props)
            end
        end

        push!(links, (gs.gene, bindings, products))
    end

    links
end

end

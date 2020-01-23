module ChainMod

using DataMod
using CellTreeMod
using ProteinPropsMod
using GeneMod
using CellMod
using ChainGraphMod

function find_chains(data::Data, ea_step::Int64, index::Int64, cell_index::Int64)
    graph = ChainGraph()
    
    DataMod.read_chunks_for_step(data, ea_step)
    for reg_step in 1:data.run.reg_steps
        tree = data.trees[ea_step][index][reg_step]
        cell = CellTreeMod.get_bf_node(tree, cell_index)
        if cell != nothing
            for gene_index in 1:length(cell.gene_states)
                gs = cell.gene_states[gene_index]
                for protein in gs.reg_site_bindings
                    if protein != nothing
                        label = ProteinPropsMod.to_str(protein.props)
                        ChainGraphMod.add_node(graph, ChainGraphMod.ProteinNode, label)
                    end
                end

                prod_rates = GeneStateMod.get_prod_rates(gs)
                for (rate, props) in zip(prod_rates, gs.prod_sites)
                    if rate != nothing
                        label = GeneMd.get_sites_str(gs.gene)
                        ChainGraphMod.add_node(graph, ChainGraphMod.GeneNode, label)
                    end
                end
            end
        end
    end
end

end

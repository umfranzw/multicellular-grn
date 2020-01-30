module ChainMod

using DataMod
using CellTreeMod
using ProteinPropsMod
using GeneMod
using CellMod
using GeneStateMod
using ChainGraphMod

function build_chain_graph(data::Data, ea_step::Int64, index::Int64, cell_index::Int64)
    graph = ChainGraph()
    
    for reg_step in 1:data.run.reg_steps
        tree = DataMod.get_tree(data, ea_step, index, reg_step)
        cell = CellTreeMod.get_bf_node(tree, cell_index)
        if cell != nothing
            for gene_index in 1:length(cell.gene_states)
                gs = cell.gene_states[gene_index]
                sites_str = GeneMod.get_sites_str(gs.gene)
                gene_label = "$(sites_str)\n($(gs.gene.genome_index))"

                bound_proteins = Set{Tuple{String, Bool}}() #(label, is_initial)
                for protein in gs.reg_site_bindings
                    if protein != nothing
                        label = ProteinPropsMod.to_str(protein.props)
                        push!(bound_proteins, (label, protein.is_initial))
                    end
                end

                prod_proteins = Set{String}() #note: produced proteins cannot be initial (so no need to track that here)
                prod_rates = GeneStateMod.get_prod_rates(gs)
                for (rate, props) in zip(prod_rates, gs.gene.prod_sites)
                    if rate != nothing
                        label = ProteinPropsMod.to_str(props)
                        push!(prod_proteins, label)
                    end
                end
                
                if !isempty(bound_proteins) || !isempty(prod_proteins)
                    #add gene node
                    ChainGraphMod.add_node(graph, ChainGraphMod.GeneNode, gene_label)
                end

                for (protein_label, is_initial) in bound_proteins
                    #add protein node
                    ChainGraphMod.add_node(graph, ChainGraphMod.ProteinNode, protein_label, is_initial)

                    #add incoming protein edge
                    ChainGraphMod.add_edge(graph, protein_label, gene_label, string(reg_step))
                end

                for protein_label in prod_proteins
                    #add protein node
                    ChainGraphMod.add_node(graph, ChainGraphMod.ProteinNode, protein_label)

                    #add outgoing protein edge
                    ChainGraphMod.add_edge(graph, gene_label, protein_label, string(reg_step))
                end
            end
        end
    end

    graph
end

end

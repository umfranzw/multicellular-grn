module ChainMod

using DataMod
using CellTreeMod
using ProteinPropsMod
using GeneMod
using CellMod
using ChainGraphMod

function build_chain_graph(data::Data, ea_step::Int64, index::Int64, cell_index::Int64)
    graph = ChainGraph()
    
    DataMod.read_chunks_for_step(data, ea_step)
    for reg_step in 1:data.run.reg_steps
        tree = data.trees[ea_step][index][reg_step]
        cell = CellTreeMod.get_bf_node(tree, cell_index)
        if cell != nothing
            for gene_index in 1:length(cell.gene_states)
                gs = cell.gene_states[gene_index]
                bound_proteins = Set{String}()
                for protein in gs.reg_site_bindings
                    if protein != nothing
                        label = ProteinPropsMod.to_str(protein.props)
                        push!(bound_proteins, label)
                        #ChainGraphMod.add_node(graph, ChainGraphMod.ProteinNode, label)
                    end
                end

                prod_proteins = Set{String}()
                prod_rates = GeneStateMod.get_prod_rates(gs)
                for (rate, props) in zip(prod_rates, gs.prod_sites)
                    if rate != nothing
                        label = ProteinPropsMod.to_str(props)
                        push!(prod_proteins, label)
                    end
                end

                sites_str = GeneMod.get_sites_str(gs.gene)
                gene_label = "$(sites_str)($(gs.gene.genome_index))"
                if !isempty(bound_proteins) || !isempty(prod_proteins)
                    #add gene node
                    ChainGraphMod.add_node(graph, ChainGraphMod.GeneNode, gene_label)
                end

                for protein_label in bound_proteins
                    #add protein node
                    ChainGraphMod.add_node(graph, ChaingGraphMod.ProteinNode, protein_label)

                    #add incoming protein edge
                    ChainGraphMod.add_edge(graph, protein_label, gene_label)
                end

                for protein_label in prod_proteins
                    #add protein node
                    ChainGraphMod.add_node(graph, ChaingGraphMod.ProteinNode, protein_label)

                    #add outgoing protein edge
                    ChainGraphMod.add_edge(graph, gene_label, protein_label)
                end
            end
        end
    end

    graph
end

end

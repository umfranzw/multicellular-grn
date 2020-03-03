module ChainGraphMod

using LightGraphs
using GraphVizMod

export ChainGraph, NodeType

@enum NodeType::UInt8 ProteinNode GeneNode

struct NodeInfo{L}
    label::L
    type::NodeType
    id::Int64
    is_initial_protein::Bool
end

#note: may add more to this later...
struct EdgeInfo{L}
    label::Union{L, Nothing}
end

mutable struct ChainGraph{L}
    graph::DiGraph
    node_info::Dict{String, NodeInfo{L}}
    edge_info::Dict{Tuple{Int64, Int64}, EdgeInfo{L}}

    function ChainGraph{L}() where L
        graph = DiGraph()
        
        new(graph, Dict{String, NodeInfo{L}}(), Dict{Tuple{Int64, Int64}, EdgeInfo{L}}())
    end
end

function add_node(graph::ChainGraph, type::NodeType, label::String, is_initial_protein::Bool=false)
    if label ∉ keys(graph.node_info)
        add_vertex!(graph.graph)
        last_id = LightGraphs.nv(graph.graph)
        graph.node_info[label] = NodeInfo(label, type, last_id, is_initial_protein)
    end
end

#note - adding the same edge twice is fine - LightGraphs will handle it
function add_edge(graph::ChainGraph, src_label::String, dest_label::String, edge_label::Union{String, Nothing}=nothing)
    src_id = graph.node_info[src_label].id
    dest_id = graph.node_info[dest_label].id
    add_edge!(graph.graph, src_id, dest_id)
    key = (src_id, dest_id)
    if key ∉ keys(graph.edge_info)
        graph.edge_info[key] = EdgeInfo(edge_label)
    end
end

function gen_dot_code(graph::ChainGraph)
    graph_buf = IOBuffer()
    print(graph_buf, "digraph G {\n")

    #generate code for nodes
    protein_buf = IOBuffer()
    gene_buf = IOBuffer()
    print(protein_buf, "{rank = same; ")
    print(gene_buf, "{rank = same; ")

    for (label, info) in graph.node_info
        if info.type == ProteinNode
            if info.is_initial_protein
                colour = "#FF0000"
            else
                colour = "#000000"
            end
            print(protein_buf, "$(info.id) [label=\"$(label)\",style=filled,fillcolor=\"#309FFF\",penwidth=4,shape=circle,color=\"$(colour)\"];\n")
            
        elseif info.type == GeneNode
            print(gene_buf, "$(info.id) [label=\"$(label)\",style=filled,fillcolor=\"#00CD66\",penwidth=4,shape=box];\n")
        end
    end
    print(protein_buf, "}\n")
    print(gene_buf, "}\n")
    
    #generate code for edges
    edge_buf = IOBuffer()
    for edge in edges(graph.graph)
        info = graph.edge_info[(edge.src, edge.dst)]
        print(edge_buf, "$(edge.src) -> $(edge.dst) [label=\"$(info.label)\"];\n")
    end

    #put it all together
    print(graph_buf, String(take!(protein_buf)))
    print(graph_buf, String(take!(gene_buf)))
    print(graph_buf, String(take!(edge_buf)))
    print(graph_buf, "}")

    String(take!(graph_buf))
end

function plot(graph::ChainGraph)
    dot_code = gen_dot_code(graph)
    
    GraphVizMod.plot(dot_code)
end

function append_for_tree(graph::ChainGraph, tree::CellTree)
    graph = ChainGraph()
    CellTreeMod.traverse(cell -> append_for_cell(graph, cell), tree.root)
end

function append_for_cell(graph::ChainGraph, cell::Cell)
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

module ChainGraphMod

using LightGraphs
using MetaGraphs
using GraphVizMod

@enum NodeType::UInt8 ProteinNode GeneNod

mutable struct ChainGraph
    graph::MetaGraph
    label_to_id::Dict{String, Int64}
end

function add_node(graph::ChainGraph, type::NodeType, label::String)
    add_vertex!(graph.graph)
    last_id = MetaGraphs.nv(graph.graph)
    set_prop!(graph.graph, last_id, :type, type)
    set_prop!(graph.graph, last_id, :label, label)
end

function add_edge(graph::ChainGraph, src::String, dest::String)
    src_id = graph.label_to_id[src]
    dest_id = graph.label_to_id[dest]
    add_edge!(graph.graph, src_id, dest_id)
end

function get_types(graph::ChainGraph)
    protein_nodes = Array{Node, 1}()
    gene_nodes = Array{Node, 1}()

    for node_id in vertices(graph.graph)
        type = get_prop(graph.graph, node_id, :type)
        if type == ProteinNode
            push!(protein_nodes, node_id)
        elseif type == GeneNode
            push!(gene_nodes, node_id)
        end
    end

    (protein_nodes, gene_nodes)
end

function gen_dot_code(graph::ChainGraph)
    protein_ids, gene_ids = get_types(graph)
    buf = IOBuffer()
    print(buf, "digraph G {\n")

    #nodes
    for ids in (protein_ids, gene_ids)
        print(buf, "{rank = same; ")
        for id in ids
            label = get_prop(graph.graph, id, :label)
            print(buf, "$(id) [label=\"$(label)\"]; ")
        end
        print(buf, "}\n")
    end

    #edges
    for edge in edges(graph.graph)
        print(buf, "$(edge.src) -> $(edge.dst);\n")
    end
    
    print(buf, "}\n")
    
    String(take!(buf))
end

function plot(graph::ChainGraph, filename::String)
    dot_code = gen_dot_code(graph)
    GraphVizMod.gen_graph(dot_code, filename)
end

end

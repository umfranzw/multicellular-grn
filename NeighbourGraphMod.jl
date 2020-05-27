module NeighbourGraphMod

using CellMod
using CellTreeMod
using ProteinPropsMod
using InteractionGraphMod
using GraphVizMod
using LightGraphs

export NeighbourGraph

mutable struct NeighbourGraph
    graph::InteractionGraph{Cell}
    tree::CellTree
    initial_props::Set{ProteinProps}

    function NeighbourGraph(tree::CellTree)
        graph = InteractionGraph{Cell}()
        CellTreeMod.traverse(cell -> InteractionGraphMod.add_obj(graph, cell), tree)
        
        new(
            graph,
            tree,
            Set{ProteinProps}()
        )
    end
end

function add_edge(graph::NeighbourGraph, src::Cell, dest::Cell, edge_label::Union{String, Nothing}=nothing)
    InteractionGraphMod.add_edge(graph.graph, src, dest, edge_label)
end

function gen_dot_code(graph::NeighbourGraph)
    graph_buf = IOBuffer()
    edge_buf = IOBuffer()
    print(graph_buf, "digraph G {\n")
    
    info = TreeInfo(graph.tree)
    for level in 1:length(info.level_to_cell)
        print(graph_buf, "{rank = same; rankdir = LR;\n")
        last_id = nothing
        for col in 1:length(info.level_to_cell[level])
            cell = info.level_to_cell[level][col]
            id = graph.graph.obj_to_id[cell]
            print(graph_buf, "$(id) [style=\"filled\",fillcolor=\"#309FFF\",penwidth=4,shape=circle,color=\"#000000\",label=\"\"];\n")

            #add invisible weighted edge between siblings
            if 2 <= col < length(info.level_to_cell[level])
                print(edge_buf, "$(last_id) -> $(id) [style=\"invis\",weight=1];\n")
            end
            last_id = id
        end
        print(graph_buf, "}\n")
    end

    #add parent-child edges
    gen_parent_child_edges(graph, graph.tree.root, edge_buf)

    #add protein edges
    for edge in edges(graph.graph.graph)
        print(edge_buf, "$(edge.src) -> $(edge.dst)")
        key = (edge.src, edge.dst)
        if key in keys(graph.graph.edge_labels)
            label = graph.graph.edge_labels[key]
            print(edge_buf, " [label=\"$(label)\"]")
        end
        print(edge_buf, ";\n")
    end

    #combine the bufs
    print(graph_buf, String(take!(edge_buf)))
    print(graph_buf, "}")

    String(take!(graph_buf))
end

function gen_parent_child_edges(graph::NeighbourGraph, cell::Cell, buf::IOBuffer)
    parent_id = graph.graph.obj_to_id[cell]
    for child in cell.children
        child_id = graph.graph.obj_to_id[child]
        print(buf, "$(parent_id) -> $(child_id);\n")
        gen_parent_child_edges(graph, child, buf)
    end
end

function plot(graph::NeighbourGraph)
    png_data = nothing
    dot_code = gen_dot_code(graph)
    png_data = GraphVizMod.plot(dot_code)

    png_data
end

end

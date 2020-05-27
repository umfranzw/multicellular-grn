module InternalGraphMod

using LightGraphs
using GraphVizMod
using GeneMod
using ProteinPropsMod
using CellTreeMod
using CellMod
using InteractionGraphMod

export InternalGraph

mutable struct InternalGraph
    graph::InteractionGraph{Union{Gene, ProteinProps}}
    initial_props::Set{ProteinProps}

    function InternalGraph()
        new(
            InteractionGraph{Union{Gene, ProteinProps}}(),
            Set{ProteinProps}()
        )
    end
end

function add_gene(graph::InternalGraph, gene::Gene)
    InteractionGraphMod.add_obj(graph.graph, gene)
end

function add_props(graph::InternalGraph, props::ProteinProps, is_initial::Bool)
    InteractionGraphMod.add_obj(graph.graph, props)
    if is_initial
        push!(graph.initial_props, props)
    end
end

function add_edge(graph::InternalGraph, src::Union{Gene, ProteinProps}, dest::Union{Gene, ProteinProps}, edge_label::Union{String, Nothing}=nothing)
    InteractionGraphMod.add_edge(graph.graph, src, dest, edge_label)
end

function gen_dot_code(graph::InternalGraph)
    graph_buf = IOBuffer()
    print(graph_buf, "digraph G {\n")

    #generate code for nodes
    protein_buf = IOBuffer()
    gene_buf = IOBuffer()
    print(protein_buf, "{rank = same; ")
    print(gene_buf, "{rank = same; rankdir = LR; ")

    genes = Array{Gene, 1}() #keep these for later (see below)
    for obj in keys(graph.graph.obj_to_id)
        id = graph.graph.obj_to_id[obj]
        if obj isa ProteinProps
            is_initial = obj in graph.initial_props
            label = ProteinPropsMod.to_str(obj, is_initial)
            
            print(protein_buf, "$(id) [label=\"$(label)\",style=filled,fillcolor=\"#309FFF\",penwidth=4,shape=circle,color=\"#000000\"];\n")
            
        elseif obj isa Gene
            push!(genes, obj)
            sites_str = GeneMod.get_sites_str(obj)
            logic = string(obj.bind_logic)
            label = "$(sites_str)\n$(logic)\n($(obj.genome_index))"
            
            print(gene_buf, "$(id) [label=\"$(label)\",style=filled,fillcolor=\"#00CD66\",penwidth=4,shape=box];\n")
        end
    end
    print(protein_buf, "}\n")
    print(gene_buf, "}\n")
    
    #generate code for edges
    edge_buf = IOBuffer()

    #draw invisible weighted edges between the genes
    #graphviz will try to keep nodes connected by weighted edges in left-to-right order if possible
    sort!(genes, by=g -> g.genome_index)
    for i in 1:length(genes) - 1
        src_id = graph.graph.obj_to_id[genes[i]]
        dest_id = graph.graph.obj_to_id[genes[i + 1]]
        print(edge_buf, "$(src_id) -> $(dest_id) [style=\"invis\",weight=1];\n")
    end
    
    for edge in edges(graph.graph.graph)
        print(edge_buf, "$(edge.src) -> $(edge.dst)")
        key = (edge.src, edge.dst)
        if key in keys(graph.graph.edge_labels)
            label = graph.graph.edge_labels[key]
            print(edge_buf, " [label=\"$(label)\"]")
        end
        print(edge_buf, ";\n")
    end

    #put it all together
    print(graph_buf, String(take!(protein_buf)))
    print(graph_buf, String(take!(gene_buf)))
    print(graph_buf, String(take!(edge_buf)))
    print(graph_buf, "}")

    String(take!(graph_buf))
end

function plot(graph::InternalGraph)
    png_data = nothing
    dot_code = gen_dot_code(graph)
    png_data = GraphVizMod.plot(dot_code)

    png_data
end

end

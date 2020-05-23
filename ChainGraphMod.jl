module ChainGraphMod

using LightGraphs
using GraphVizMod
using GeneMod
using ProteinMod
using ProteinPropsMod
using CellTreeMod
using CellMod
using DataStructures

export ChainGraph

mutable struct ChainGraph
    graph::DiGraph
    id_to_obj::Dict{Int64, Union{Gene, ProteinProps}}
    initial_props::Set{ProteinProps}
    obj_to_id::OrderedDict{Union{Gene, ProteinProps}, Int64}
    #note: it's possible for unlabelled edges to exist.
    #In that case, the edge will be in graph, but not in this dict.
    edge_labels::Dict{Tuple{Int64, Int64}, String}

    function ChainGraph()
        graph = DiGraph()
        
        new(
            graph,
            Dict{Int64, Union{Gene, ProteinProps}}(),
            Set{ProteinProps}(),
            OrderedDict{Union{Gene, ProteinProps}, Int64}(),
            Dict{Tuple{Int64, Int64}, String}()
        )
    end
end

function add_gene(graph::ChainGraph, gene::Gene)
    add_node(graph, gene)
end

function add_props(graph::ChainGraph, props::ProteinProps, is_initial::Bool)
    add_node(graph, props)
    if is_initial
        push!(graph.initial_props, props)
    end
end

#don't call this one - call one of the above two
function add_node(graph::ChainGraph, node::Union{Gene, ProteinProps})
    if node ∉ keys(graph.obj_to_id)
        add_vertex!(graph.graph)
        last_id = LightGraphs.nv(graph.graph)
        graph.id_to_obj[last_id] = node
        graph.obj_to_id[node] = last_id
    end
end

function add_edge(graph::ChainGraph, src::Union{Gene, ProteinProps}, dest::Union{Gene, ProteinProps}, edge_label::Union{String, Nothing}=nothing)
    src_id = graph.obj_to_id[src]
    dest_id = graph.obj_to_id[dest]

    #note - adding the same edge twice is fine - LightGraphs will handle it
    add_edge!(graph.graph, src_id, dest_id)
    if edge_label != nothing
        graph.edge_labels[(src_id, dest_id)] = edge_label
    end
end

function gen_dot_code(graph::ChainGraph)
    graph_buf = IOBuffer()
    print(graph_buf, "digraph G {\n")

    #generate code for nodes
    protein_buf = IOBuffer()
    gene_buf = IOBuffer()
    print(protein_buf, "{rank = same; ")
    print(gene_buf, "{rank = same; rankdir = LR; ")

    genes = Array{Gene, 1}() #keep these for later (see below)
    for obj in keys(graph.obj_to_id)
        id = graph.obj_to_id[obj]
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
        src_id = graph.obj_to_id[genes[i]]
        dest_id = graph.obj_to_id[genes[i + 1]]
        print(edge_buf, "$(src_id) -> $(dest_id) [style=\"invis\",weight=1];\n")
    end
    
    for edge in edges(graph.graph)
        print(edge_buf, "$(edge.src) -> $(edge.dst)")
        key = (edge.src, edge.dst)
        if key in keys(graph.edge_labels)
            label = graph.edge_labels[key]
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

function plot(graph::ChainGraph)
    png_data = nothing
    dot_code = gen_dot_code(graph)
    png_data = GraphVizMod.plot(dot_code)

    png_data
end

end

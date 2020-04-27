module ChainGraphMod

using LightGraphs
using GraphVizMod
using GeneMod
using ProteinMod
using ProteinPropsMod
using CellTreeMod
using CellMod

export ChainGraph

mutable struct ChainGraph
    graph::DiGraph
    id_to_obj::Dict{Int64, Union{Gene, Tuple{ProteinProps, Bool}}}
    obj_to_id::Dict{Union{Gene, Tuple{ProteinProps, Bool}}, Int64}
    #note: it's possible for unlabelled edges to exist.
    #In that case, the edge will be in graph, but not in this dict.
    edge_labels::Dict{Tuple{Int64, Int64}, String}

    function ChainGraph()
        graph = DiGraph()
        
        new(
            graph,
            Dict{Int64, Union{Gene, Tuple{ProteinProps, Bool}}}(),
            Dict{Union{Gene, Tuple{ProteinProps, Bool}}, Int64}(),
            Dict{Tuple{Int64, Int64}, String}()
        )
    end
end

function is_empty(graph::ChainGraph)
    length(graph.id_to_obj) == 0
end

function add_node(graph::ChainGraph, node::Union{Gene, Tuple{ProteinProps, Bool}})
    if node ∉ keys(graph.obj_to_id)
        add_vertex!(graph.graph)
        last_id = LightGraphs.nv(graph.graph)
        graph.id_to_obj[last_id] = node
        graph.obj_to_id[node] = last_id
    end
end

function add_edge(graph::ChainGraph, src::Union{Gene, Tuple{ProteinProps, Bool}}, dest::Union{Gene, Tuple{ProteinProps, Bool}}, edge_label::Union{String, Nothing}=nothing)
    src_id = graph.obj_to_id[src]
    dest_id = graph.obj_to_id[dest]

    #note - adding the same edge twice is fine - LightGraphs will handle it
    add_edge!(graph.graph, src_id, dest_id)
    if edge_label != nothing
        graph.edge_labels[(src_id, dest_id)] = edge_label
    end
end

# function get_app_contributing_genes(graph::ChainGraph)
#     #Work backwards, starting from the app proteins and going up to the initial_proteins.
#     #Record all genes involved along the way.
#     #Don't forget to account for cycles!
#     app_contrib_genes = Set{Gene}()

#     #note: we start and end with protein ids
#     #grab all of the app protein ids
#     protein_ids = Set{Int64}()
#     for (id, obj) in graph.id_to_obj
#         if obj isa Protein && obj.props.type == ProteinPropsMod.Application
#             push!(protein_ids, id)
#         end
#     end

#     #go back one step at a time (from protein -> gene -> protein)
#     while !isempty(protein_ids)
#         gene_ids = Set{Int64}()
#         for protein_id in protein_ids
#             #get ids of all (gene) nodes with an edge to node with protein_id
#             src_ids = LightGraphs.inneighbors(protein_id)

#             #avoid cycles - filter out all app-contributing genes that we've already encountered
#             for id in src_ids
#                 gene = graph.id_to_obj[id]
#                 if gene ∉ app_contrib_genes
#                     push!(gene_ids, id)
#                     push!(app_contrib_genes, gene)
#                 end
#             end
#         end

#         protein_ids = Set{Int64}()
#         for gene_id in gene_ids
#             #get ids of all (protein) nodes with an edge to node with gene_id
#             src_ids = LightGraph.inneighbors(gene_id)

#             #filter out all the initial proteins (since there's nothing before them)
#             for id in src_ids
#                 if !graph.id_to_obj[id].is_initial
#                     push!(protein_ids, id)
#                 end
#             end
#         end
#     end

#     app_contrib_genes
# end

function gen_dot_code(graph::ChainGraph)
    graph_buf = IOBuffer()
    print(graph_buf, "digraph G {\n")

    #generate code for nodes
    protein_buf = IOBuffer()
    gene_buf = IOBuffer()
    print(protein_buf, "{rank = same; ")
    print(gene_buf, "{rank = same; ")

    for obj in keys(graph.obj_to_id)
        id = graph.obj_to_id[obj]
        if obj isa Tuple
            props, is_initial = obj
            label = ProteinPropsMod.to_str(props)
            
            if is_initial
                colour = "#FF0000"
            else
                colour = "#000000"
            end
            print(protein_buf, "$(id) [label=\"$(label)\",style=filled,fillcolor=\"#309FFF\",penwidth=4,shape=circle,color=\"$(colour)\"];\n")
            
        elseif obj isa Gene
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
    if !is_empty(graph)
        dot_code = gen_dot_code(graph)
        png_data = GraphVizMod.plot(dot_code)
    end

    png_data
end

# function append_for_tree(graph::ChainGraph, tree::CellTree, reg_step::Int64)
#     graph = ChainGraph()
#     CellTreeMod.traverse(cell -> append_for_cell(graph, cell, reg_step), tree.root)
# end

# function append_for_cell(graph::ChainGraph, cell::Cell, reg_step::Int64)
#     for gene_index in 1:length(cell.gene_states)
#         gs = cell.gene_states[gene_index]

#         bound_proteins = Set{Protein}()
#         for protein in gs.reg_site_bindings
#             if protein != nothing
#                 push!(bound_proteins, protein)
#             end
#         end

#         prod_proteins = Set{Protein}() #note: produced proteins cannot be initial (so no need to track that here)
#         prod_rates = GeneStateMod.get_prod_rates(gs)
#         for (rate, props) in zip(prod_rates, gs.gene.prod_sites)
#             if rate != nothing
#                 protein = ProteinStoreMod.get(cell.proteins, props)
#                 push!(prod_proteins, protein)
#             end
#         end
        
#         if !isempty(bound_proteins) || !isempty(prod_proteins)
#             #add gene node
#             ChainGraphMod.add_node(graph, ChainGraphMod.GeneNode, gs.gene)
#         end

#         for protein in bound_proteins
#             #add protein node
#             ChainGraphMod.add_node(graph, ChainGraphMod.ProteinNode, (protein.props, protein.is_initial))

#             #add incoming protein edge
#             ChainGraphMod.add_edge(graph, protein, gs.gene, string(reg_step))
#         end

#         for protein in prod_proteins
#             #add protein node
#             ChainGraphMod.add_node(graph, ChainGraphMod.ProteinNode, (protein.props, protein.is_initial))

#             #add outgoing protein edge
#             ChainGraphMod.add_edge(graph, gs.gene, (protein.props, protein.is_initial), string(reg_step))
#         end
#     end
# end

end

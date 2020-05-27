module InteractionGraphMod

using LightGraphs
using DataStructures

export InteractionGraph

mutable struct InteractionGraph{T}
    graph::DiGraph
    id_to_obj::Dict{Int64, T}
    obj_to_id::Dict{T, Int64}
    #note: it's possible for unlabelled edges to exist.
    #In that case, the edge will be in graph, but not in this dict.
    edge_labels::Dict{Tuple{Int64, Int64}, String}

    function InteractionGraph{T}() where T
        new(
            DiGraph(),
            Dict{Int64, T}(),
            OrderedDict{T, Int64}(),
            Dict{Tuple{Int64, Int64}, String}()
        )
    end
end

function add_obj(graph::InteractionGraph{T}, obj::T) where T
    if obj ∉ keys(graph.obj_to_id)
        add_vertex!(graph.graph)
        last_id = LightGraphs.nv(graph.graph)
        graph.id_to_obj[last_id] = obj
        graph.obj_to_id[obj] = last_id
    end
end

function add_edge(graph::InteractionGraph{T}, src::T, dest::T, edge_label::Union{String, Nothing}=nothing) where T
    src_id = graph.obj_to_id[src]
    dest_id = graph.obj_to_id[dest]

    #note - adding the same edge twice is fine - LightGraphs will handle it
    add_edge!(graph.graph, src_id, dest_id)
    if edge_label != nothing
        graph.edge_labels[(src_id, dest_id)] = edge_label
    end
end

end

module GeneStateMod

using GeneMod
using ProteinMod
using RunMod

export GeneState

mutable struct GeneState
    gene::Gene
    output_rate::Float64
    bind_site_bindings::Array{Union{Protein, Nothing}, 1}
    prod_site_bindings::Array{Union{Protein, Nothing}, 1}
    active_products::Dict{BitArray, Protein}

    function GeneState(run::Run, gene::Gene)
        bind_site_bindings = []
        new(
            gene,
            0.0,
            fill(nothing, run.num_bind_sites),
            fill(nothing, run.num_bind_sites),
            Dict{BitArray, Protein}()
        )
    end
end

function bind(gs::GeneState, protein::Protein, site_index::Int64)
    gs.bind_site_bindings[site_index] = protein

    if protein.seq ∉ gs.active_products
        gs.active_products[protein.seq] = protein
    end
end

function unbind(gs::GeneState, site_index::Int64)
    gs.bind_site_bindings[site_index] = nothing
end

end

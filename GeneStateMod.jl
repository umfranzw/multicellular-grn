module GeneStateMod

@enum SiteType::Bool BindSite=false ProdSite=true

using GeneMod
using ProteinMod
using RunMod

export GeneState, SiteType

mutable struct GeneState
    run::Run
    gene::Gene
    output_rate::Float64
    bind_site_bindings::Array{Union{Protein, Nothing}, 1}
    prod_site_bindings::Array{Union{Protein, Nothing}, 1}

    function GeneState(run::Run, gene::Gene)
        bind_site_bindings = []
        new(
            run,
            gene,
            0.0,
            fill(nothing, run.num_bind_sites),
            fill(nothing, run.num_bind_sites)
        )
    end
end

function bind(gs::GeneState, protein::Protein, site_type::SiteType, site_index::Int64)
    if site_type == BindSite
        gs.bind_site_bindings[site_index] = protein
    else
        gs.prod_site_bindings[site_index] = protein
    end
end

function unbind(gs::GeneState, site_type::SiteType, site_index::Int64)
    if site_type == BindSite
        gs.bind_site_bindings[site_index] = nothing
    else
        gs.prod_site_bindings[site_index] = nothing
    end
end

function get_active_prod_seqs(gs::GeneState)
    seqs = []
    for i in 1:gs.run.num_bind_sites
        if gs.bind_site_bindings[i] != nothing && gs.prod_site_bindings == nothing
            append!(seqs, gs.gene.prod_sites[i])
        end
    end
    
    seqs
end

end

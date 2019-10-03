module GeneStateMod

@enum SiteType::Int64 RegSite GrowthSite BindSite ProdSite

using GeneMod
using ProteinMod
using RunMod

export GeneState, SiteType

mutable struct GeneState
    run::Run
    gene::Gene
    reg_site_bindings::Array{Union{Protein, Nothing}, 1}
    prod_site_bindings::Array{Union{Protein, Nothing}, 1}

    function GeneState(run::Run, gene::Gene)
        new(
            run,
            gene,
            repeat([nothing], length(instances(GeneMod.RegSites)))
            repeat([nothing], length(instances(GeneMod.ProdSites)))
        )
    end
end

function bind(gs::GeneState, protein::Protein, site::Union{GeneMod.RegSites, GeneMod.ProdSites})
    modify_binding(gs, site, protein)
end

function unbind(gs::GeneState, site::Union{GeneMod.RegSites, GeneMod.ProdSites})
    modify_binding(gs, site, nothing)
end

function modify_binding(gs::GeneState, site::Union{GeneMod.RegSites, GeneMod.ProdSites}, val::Union{Protein, Nothing})
    index = Int64(site)
    if site isa ProteinMod.RegSites
        reg_site_bindings[index] = val
    else
        prod_site_bindings[index] = val
    end
end

function get_active_prod_props(gs::GeneState)
    props = []
    #TODO
    #combine influences of both intra sites and use them to activate the production site
    #do same for inter sites
    
    seqs
end

end

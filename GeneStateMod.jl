module GeneStateMod

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
    if site isa GeneMod.RegSites
        reg_site_bindings[index] = val
    else
        prod_site_bindings[index] = val
    end
end

function get_binding_state(gs::GeneState, site::Union{GeneMod.RegSites, GeneMod.ProdSites})
    index = Int64(site)
    if site isa GeneMod.RegSites
        return reg_site_bindings[index]
    else
        return prod_site_bindings[index]
    end
end

function calc_rate_for_sites(gs::GeneState, reg_sites::Array{GeneMod.RegSites, 1})
    weight = 0.0
    for site in reg_sites
        protein = get_binding_state(gs, site)
        if protein != nothing && protein.props.reg_action != ProteinMod.Inhibit
            weight += protein.concs[gs.gene.genome_index]
        end
    end
    
    weight / length(reg_sites) #take the average
end

function get_prod_rates(gs::GeneState)
    #calculate the influence of each pair of reg sites on the prod site they regulate
    #this is a value in [0.0, 1.0]
    
    #For intra production site
    intra_rate = calc_rate_for_sites(gs, [GeneMod.IntraIntra, GeneMod.InterIntra])
    
    #for inter production site
    inter_rate = calc_rate_for_sites(gs, [GeneMod.IntraInter, GeneMod.InterInter])

    (intra=intra_rate, inter=inter_rate)
end

end

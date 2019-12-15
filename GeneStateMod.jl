module GeneStateMod

using GeneMod
using ProteinMod
using ProteinPropsMod
using RunMod
using MiscUtilsMod
using CustomEnumMod

import Base.show

export GeneState

mutable struct GeneState
    config::Config
    gene::Gene
    reg_site_bindings::Array{Union{Protein, Nothing}, 1}
    prod_site_bindings::Array{Union{Protein, Nothing}, 1}

    function GeneState(config::Config, gene::Gene)
        new(
            config,
            gene,
            repeat([nothing], length(GeneMod.RegSites)),
            repeat([nothing], length(GeneMod.ProdSites))
        )
    end
end

function show(io::IO, gs::GeneState, ilevel::Int64=0)
    iprintln(io, "GeneState", ilevel)

    iprintln(io, "genome_index: $(gs.gene.genome_index)", ilevel + 1)
    iprintln(io, "reg_site_bindings", ilevel + 1)
    iprint(io, "", ilevel + 2)
    for (sym, val) in GeneMod.RegSites
        print(io, "$(string(sym)): ")
        
        site = gs.reg_site_bindings[val]
        if site == nothing
            print(io, "(nothing)  ")
        else
            print(io, site.props)
            print(io, "  ")
        end
    end
    println(io, "")

    iprintln(io, "prod_site_bindings", ilevel + 1)
    iprint(io, "", ilevel + 2)
    for (sym, val) in GeneMod.ProdSites
        print(io, "$(string(sym)): ")
        
        site = gs.prod_site_bindings[val]
        if site == nothing
            print(io, "(nothing)  ")
        else
            print(io, site.props)
            print(io, "  ")
        end
    end
    println(io, "")
end

function bind(gs::GeneState, protein::Protein, site::CustomVal)
    modify_binding(gs, site, protein)
end

function unbind(gs::GeneState, site::CustomVal)
    modify_binding(gs, site, nothing)
end

function modify_binding(gs::GeneState, site::CustomVal, val::Union{Protein, Nothing})
    index = Int64(site)
    if site isa CustomEnumMod.RegSiteVal
        gs.reg_site_bindings[index] = val
    else
        gs.prod_site_bindings[index] = val
    end
end

function get_binding_state(gs::GeneState, site::CustomVal)
    index = Int64(site)
    if site isa CustomEnumMod.RegSiteVal
        return gs.reg_site_bindings[index]
    else
        return gs.prod_site_bindings[index]
    end
end

#reg_sites is array of enum values from GeneMod.RegSites
function calc_rate_for_reg_sites(gs::GeneState, sites::Array{CustomEnumMod.RegSiteVal, 1})
    producing = false
    weight = 0.0
    for site in GeneMod.RegSites
        protein = get_binding_state(gs, site)
        if protein != nothing && protein.props.reg_action != ProteinPropsMod.Inhibit
            producing = true
            weight += protein.concs[gs.gene.genome_index]
        end
    end

    if producing
        return weight / length(GeneMod.RegSites) #take the average (this will be a value in [0.0, 1.0])
    else
        return nothing
    end
end

function get_prod_rates(gs::GeneState)
    #calculate the influence of each pair of reg sites on the prod site they regulate
    #this is a value in [0.0, 1.0]
    
    #For intra production site
    intra_rate = calc_rate_for_reg_sites(gs, [GeneMod.RegSites.IntraIntra, GeneMod.RegSites.InterIntra])
    
    #for inter production site
    inter_rate = calc_rate_for_reg_sites(gs, [GeneMod.RegSites.IntraInter, GeneMod.RegSites.InterInter])

    (intra=intra_rate, inter=inter_rate)
end

end

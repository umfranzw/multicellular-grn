module GeneStateMod

using GeneMod
using ProteinMod
using ProteinPropsMod
using RunMod
using MiscUtilsMod

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
            repeat([nothing], length(instances(GeneMod.RegSites))),
            repeat([nothing], length(instances(GeneMod.ProdSites)))
        )
    end
end

function show(io::IO, gs::GeneState, ilevel::Int64=0)
    iprintln(io, "GeneState", ilevel)

    iprintln(io, "genome_index: $(gs.gene.genome_index)", ilevel + 1)
    iprintln(io, "reg_site_bindings", ilevel + 1)
    iprint(io, "", ilevel + 2)
    for site_type in instances(GeneMod.RegSites)
        print(io, "$(string(site_type)): ")
        
        site = gs.reg_site_bindings[Int64(site_type)]
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
    for site_type in instances(GeneMod.ProdSites)
        print(io, "$(string(site_type)): ")
        
        site = gs.prod_site_bindings[Int64(site_type)]
        if site == nothing
            print(io, "(nothing)  ")
        else
            print(io, site.props)
            print(io, "  ")
        end
    end
    println(io, "")
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
        gs.reg_site_bindings[index] = val
    else
        gs.prod_site_bindings[index] = val
    end
end

function get_binding_state(gs::GeneState, site::Union{GeneMod.RegSites, GeneMod.ProdSites})
    index = Int64(site)
    if site isa GeneMod.RegSites
        return gs.reg_site_bindings[index]
    else
        return gs.prod_site_bindings[index]
    end
end

function calc_rate_for_sites(gs::GeneState, reg_sites::Array{GeneMod.RegSites, 1})
    producing = false
    weight = 0.0
    for site in reg_sites
        protein = get_binding_state(gs, site)
        if protein != nothing && protein.props.reg_action != ProteinPropsMod.Inhibit
            producing = true
            weight += protein.concs[gs.gene.genome_index]
        end
    end

    if producing
        return weight / length(reg_sites) #take the average (this will be a value in [0.0, 1.0])
    else
        return nothing
    end
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

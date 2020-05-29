module GeneMod

using RunMod
using ProteinPropsMod
using MiscUtilsMod
using Printf

import RandUtilsMod
import Random
import Base.show

export Gene, BindSite, ProdSite

@enum BindLogic::UInt8 Id=1 And Or Xor

mutable struct BindSite
    #from the protein
    type::ProteinPropsMod.ProteinType
    #no fcn here
    action::ProteinPropsMod.ProteinAction
    loc::ProteinPropsMod.ProteinLoc
    #additional binding params
    threshold::Float64
    consum_rate::Float64
end

mutable struct ProdSite
    type::ProteinPropsMod.ProteinType
    fcn::ProteinPropsMod.ProteinFcn
    action::ProteinPropsMod.ProteinAction
    loc::ProteinPropsMod.ProteinLoc
    #no arg here - it's a binding "parameter"
end

mutable struct Gene
    config::Config
    genome_index::Int64
    bind_logic::BindLogic
    bind_sites::Array{BindSite, 1}
    prod_sites::Array{ProdSite, 1}
end

function get_sites_str(gene::Gene)
    buf = IOBuffer()
    
    for i in 1:length(gene.bind_sites)
        bind_str = GeneMod.get_bind_site_str(gene, i)
        print(buf, bind_str)
        
        if i < length(gene.bind_sites)
            print(buf, ", ")
        end
    end
    
    print(buf, " : ")

    for i in 1:length(gene.prod_sites)
        prod_str = GeneMod.get_prod_site_str(gene, i)
        print(buf, prod_str)
        if i < length(gene.prod_sites)
            print(buf, ", ")
        end
    end
    
    String(take!(buf))
end

function get_bind_site_str(gene::Gene, index::Int64)
    site = gene.bind_sites[index]
    desc = ProteinPropsMod.get_abbrev_props_str(
        type=site.type,
        action=site.action,
        loc=site.loc
    )
    #@printf(buf, " %0.1f", gene.bind_sites[index].threshold)
    #@printf(buf, " %0.1f", gene.bind_sites[index].consum_rate)

    desc
end

function get_prod_site_str(gene::Gene, index::Int64)
    site = gene.prod_sites[index]
    desc = ProteinPropsMod.get_abbrev_props_str(
        type=site.type,
        fcn=site.fcn,
        action=site.action,
        loc=site.loc
    )

    desc
end

function show(io::IO, gene::Gene, ilevel::Int64=0)
    iprintln(io, "Gene:", ilevel)
    iprintln(io, "genome_index: $(gene.genome_index)", ilevel + 1)
    iprintln(io, "bind_logic: $(string(gene.bind_logic))", ilevel + 1)
    
    iprintln(io, "sites:", ilevel + 1)
    iprint(io, get_sites_str(gene), ilevel + 2)
end

function rand_prod_site(
    config::Config;
    type::Union{Array{ProteinPropsMod.ProteinType, 1}, Nothing}=nothing,
    fcn::Union{Array{ProteinPropsMod.ProteinFcn, 1}, Nothing}=nothing,
    action::Union{Array{ProteinPropsMod.ProteinAction, 1}, Nothing}=nothing,
    loc::Union{Array{ProteinPropsMod.ProteinLoc, 1}, Nothing}=nothing
)
    args = Array{Any, 1}()
    pairs = (
        (ProteinPropsMod.ProteinType, type),
        (ProteinPropsMod.ProteinFcn, fcn),
        (ProteinPropsMod.ProteinAction, action),
        (ProteinPropsMod.ProteinLoc, loc)
    )
    foreach(p -> push!(args, ProteinPropsMod.rand_prop(config, p...)), pairs)

    ProdSite(args...)
end

function rand_bind_site(
    config::Config;
    type::Union{Array{ProteinPropsMod.ProteinType, 1}, Nothing}=nothing,
    action::Union{Array{ProteinPropsMod.ProteinAction, 1}, Nothing}=nothing,
    loc::Union{Array{ProteinPropsMod.ProteinLoc, 1}, Nothing}=nothing,
    threshold::Union{Array{Float64, 1}, Nothing}=nothing,
    consum_rate::Union{Array{Float64, 1}, Nothing}=nothing
)
    args = Array{Any, 1}()
    pairs = (
        (ProteinPropsMod.ProteinType, type),
        (ProteinPropsMod.ProteinAction, action),
        (ProteinPropsMod.ProteinLoc, loc)
    )
    foreach(p -> push!(args, ProteinPropsMod.rand_prop(config, p...)), pairs)
    
    if threshold == nothing
        threshold_val = RandUtilsMod.rand_float(config)
    else
        threshold_val = Random.rand(config.rng, threshold)
    end
    push!(args, threshold_val)
    
    if consum_rate == nothing
        consum_rate_val = RandUtilsMod.rand_float(config)
    else
        consum_rate_val = Random.rand(config.rng, consum_rate)
    end
    push!(args, consum_rate_val)

    BindSite(args...)
end

#note: binding sites should always be of type Internal
function rand_init(
    config::Config,
    genome_index::Int64,
    bind_site_types::Array{ProteinPropsMod.ProteinType, 1}=[ProteinPropsMod.Internal],
    prod_site_types::Array{ProteinPropsMod.ProteinType, 1}=[],
    bind_logic::Union{Array{BindLogic, 1}, Nothing}=nothing
)
    bind_sites = Array{BindSite, 1}()
    prod_sites = Array{ProdSite, 1}()
    for i in 1:config.run.bind_sites_per_gene
        bind_site = rand_bind_site(
            config,
            type=bind_site_types
        )
        push!(bind_sites, bind_site)

        prod_site = rand_prod_site(
            config
        )
        push!(prod_sites, prod_site)
    end

    if bind_logic == nothing
        logic = RandUtilsMod.rand_enum_val(config, BindLogic)
    else
        logic = Random.rand(config.rng, bind_logic)
    end
    
    Gene(
        config,
        genome_index,
        logic,
        bind_sites,
        prod_sites
    )
end

end

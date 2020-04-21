module GeneMod

using RunMod
using ProteinPropsMod
using MiscUtilsMod
using Printf

import RandUtilsMod
import Random
import Base.show

export Gene, BindSite

@enum BindLogic::Int8 Id=1 And Or Xor

mutable struct BindSite
    #from the protein
    type::ProteinPropsMod.ProteinType
    action::ProteinPropsMod.ProteinAction
    loc::ProteinPropsMod.ProteinLoc
    #additional binding params
    threshold::Float64
    consum_rate::Float64
end

mutable struct Gene
    config::Config
    genome_index::Int64
    bind_logic::BindLogic
    bind_sites::Array{BindSite, 1}
    prod_sites::Array{ProteinProps, 1}
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
        prod_str = GeneMode.get_prod_site_str(gene, i)
        print(buf, prod_str)
        if i < length(gene.prod_sites)
            print(buf, ", ")
        end
    end
    
    String(take!(buf))
end

function get_bind_site_str(gene::Gene, index::Int64)
    buf = IOBuffer()
    pairs = (
        (ProteinPropsMod.ProteinType, gene.bind_sites[index].type),
        (ProteinPropsMod.ProteinAction, gene.bind_sites[index].action),
        (ProteinPropsMod.ProteinLoc, gene.bind_sites[index].loc)
    )
    for (enum, val) in pairs
        chunk = string(val)[1:3]
        print(buf, chunk)
    end
    #@printf(buf, " %0.1f", gene.bind_sites[index].threshold)
    #@printf(buf, " %0.1f", gene.bind_sites[index].consum_rate)

    String(take!(buf))
end

function get_prod_site_str(gene::Gene, index::Int64)
    ProteinPropsMod.to_str(gene.prod_sites[index])
end

function show(io::IO, gene::Gene, ilevel::Int64=0)
    iprintln(io, "Gene:", ilevel)
    iprintln(io, "genome_index: $(gene.genome_index)", ilevel + 1)
    iprintln(io, "bind_logic: $(string(gene.bind_logic))", ilevel + 1)
    
    iprintln(io, "sites:", ilevel + 1)
    iprint(io, get_sites_str(gene), ilevel + 2)
end

function rand_bind_site(
    config::Config;
    type::Union{Array{ProteinPropsMod.ProteinType, 1}, Nothing}=nothing,
    loc::Union{Array{ProteinPropsMod.ProteinLoc, 1}, Nothing}=nothing,
    action::Union{Array{ProteinPropsMod.ProteinAction, 1}, Nothing}=nothing,
    threshold::Union{Array{Float64, 1}, Nothing}=nothing,
    consum_rate::Union{Array{Float64, 1}, Nothing}=nothing
)
    if type == nothing
        type_val = RandUtilsMod.rand_enum_val(config, ProteinPropsMod.ProteinType)
    else
        type_val = Random.rand(config.rng, type)
    end
    if action == nothing
        action_val = RandUtilsMod.rand_enum_val(config, ProteinPropsMod.ProteinAction)
    else
        action_val = Random.rand(config.rng, type)
    end
    if loc == nothing
        loc_val = RandUtilsMod.rand_enum_val(config, ProteinPropsMod.ProteinLoc)
    else
        loc_val = Random.rand(config.rng, loc)
    end
    if threshold == nothing
        threshold_val = RandUtilsMod.rand_float(config)
    else
        threshold_val = Random.rand(config.rng, threshold)
    end
    if consum_rate == nothing
        consum_rate_val = RandUtilsMod.rand_float(config)
    else
        consum_rate_val = Random.rand(config.rng, consum_rate)
    end

    BindSite(type_val, action_val, loc_val, threshold_val, consum_rate_val)
end

function rand_init(config::Config, genome_index::Int64)
    bind_sites = Array{BindSite, 1}()
    prod_sites = Array{ProteinProps, 1}()
    for i in 1:config.run.bind_sites_per_gene
        bind_site = rand_bind_site(
            config,
            #note: binding sites should never be of type Application
            type=[ProteinPropsMod.Internal, ProteinPropsMod.Neighbour, ProteinPropsMod.Diffusion]
        )
        push!(bind_sites, bind_site)

        prod_site = ProteinPropsMod.rand_init(
            config
        )
        push!(prod_sites, prod_site)
    end

    bind_logic = RandUtilsMod.rand_enum_val(config, BindLogic)
    
    Gene(
        config,
        genome_index,
        bind_logic,
        bind_sites,
        prod_sites
    )
end

end

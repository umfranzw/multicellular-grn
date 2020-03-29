module GeneMod

using RunMod
using ProteinPropsMod
using MiscUtilsMod

import RandUtilsMod
import Random
import Base.show
import Formatting

export Gene

@enum BindLogic::Int8 Id=1 And Or Xor

mutable struct BindSite
    #from the protein
    type::ProteinProps.ProteinType
    action::ProteinProps.ProteinAction
    loc::ProteinProps.ProteinLoc
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
        pairs = (
            (ProteinProps.ProteinType, gene.bind_sites[i].type),
            (ProteinProps.ProteinAction, gene.bind_sites[i].action),
            (ProteinProps.ProteinLoc, gene.bind_sites[i].loc)
        )
        for (enum, val) in pairs
            width = MiscUtilsMod.digits_needed(length(instances(enum)))
            fs = Formatting.FormatSpec("0$(width)d")
            str = Formatting.fmt(fs, Int64(val))
            print(buf, str)
        end
        @printf(buf, " %0.1f", gene.binding_sites[i].threshold)
        @printf(buf, " %0.1f", gene.binding_sites[i].consum_rate)
        
        if i < length(gene.reg_sites)
            print(buf, ", ")
        end
    end
    print(buf, "\n")

    for i in 1:length(gene.prod_sites)
        print(buf, ProteinPropsMod.to_str(gene.prod_sites[i]))
        if i < length(gene.prod_sites)
            print(buf, ", ")
        end
    end
    
    String(take!(buf))
end

function show(io::IO, gene::Gene, ilevel::Int64=0)
    iprintln(io, "Gene:", ilevel)
    iprintln(io, "genome_index: $(gene.genome_index)", ilevel + 1)
    
    iprintln(io, "sites:", ilevel + 1)
    iprint(io, get_sites_str(gene), ilevel + 2)
end
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
        threshold_val = RandUtilsMod.rand_float(config.rng)
    else
        threshold_val = Random.rand(config.rng, threshold)
    end
    if consum_rate == nothing
        consum_rate_val = RandUtilsMod.rand_float(config.rng)
    else
        consum_rate_val = Random.rand(config.rng, consum_rate)
    end

    BindSite(type_val, loc_val, threshold_val, consum_rate_val)
end

function rand_init(config::Config, genome_index::Int64)
    bind_sites = Array{BindSite, 1}()
    prod_sites = Array{ProteinProps, 1}()
    for i in 1:config.run.num_bind_sites
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

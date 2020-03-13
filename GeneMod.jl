module GeneMod

using RunMod
using ProteinPropsMod
using MiscUtilsMod

import RandUtilsMod
import Random
import Base.show

export Gene

@enum RegSite::Int8 IntraIntra=1 IntraInter InterIntra InterInter

@enum ProdSite::Int8 Intra=1 Inter

mutable struct Gene
    config::Config
    genome_index::Int64
    reg_sites::Array{ProteinProps, 1}
    prod_sites::Array{ProteinProps, 1}
end

function get_sites_str(gene::Gene)
    buf = IOBuffer()
    for i in 1:length(gene.reg_sites)
        print(buf, ProteinPropsMod.to_str(gene.reg_sites[i]))
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
    
    iprintln(io, "reg_sites:", ilevel + 1)
    for val in instances(RegSite)
        iprint(io, "$(string(val)): ", ilevel + 2)
        ProteinPropsMod.show(io, gene.reg_sites[Int64(val)], 0)
    end

    iprintln(io, "prod_sites:", ilevel + 1)
    for val in instances(ProdSite)
        iprint(io, "$(string(val)): ", ilevel + 2)
        ProteinPropsMod.show(io, gene.prod_sites[Int64(val)], 0)
    end
end

function rand_site(
    config::Config;
    type::Union{Array{ProteinPropsMod.ProteinType, 1}, Nothing}=nothing,
    target::Union{Array{ProteinPropsMod.ProteinTarget, 1}, Nothing}=nothing,
    reg_action::Union{Array{ProteinPropsMod.ProteinRegAction, 1}, Nothing}=nothing,
    app_action::Union{Array{UInt8, 1}, Nothing}=nothing
)
    if type == nothing
        type_val = RandUtilsMod.rand_enum_val(config, ProteinPropsMod.ProteinType)
    else
        type_val = Random.rand(config.rng, type)
    end
    if target == nothing
        target_val = RandUtilsMod.rand_enum_val(config, ProteinPropsMod.ProteinTarget)
    else
        target_val = Random.rand(config.rng, target)
    end
    if reg_action == nothing
        reg_action_val = RandUtilsMod.rand_enum_val(config, ProteinPropsMod.ProteinRegAction)
    else
        reg_action_val = Random.rand(config.rng, reg_action)
    end
    if app_action == nothing
        app_action_val = UInt8(RandUtilsMod.rand_int(config, 1, Int64(ProteinPropsMod.num_app_actions)))
    else
        app_action_val = Random.rand(config.rng, app_action)
    end

    ProteinProps(type_val, target_val, reg_action_val, app_action_val)
end

function rand_init(config::Config, genome_index::Int64)
    #reg_sites
    reg_sites = Array{ProteinProps, 1}()
    intra_intra = rand_site(
        config,
        type=ProteinPropsMod.Reg,
        target=ProteinPropsMod.Intra
    )
    push!(reg_sites, intra_intra)

    intra_inter = rand_site(
        config,
        type=ProteinPropsMod.Reg,
        target=ProteinPropsMod.Intra
    )
    push!(reg_sites, intra_inter)

    inter_intra = rand_site(
        config,
        type=ProteinPropsMod.Reg,
        target=[ProteinPropsMod.InterLocal, ProteinPropsMod.InterDistant]
    )
    push!(reg_sites, inter_intra)

    inter_inter = rand_site(
        config,
        type=ProteinPropsMod.Reg,
        target=[ProteinPropsMod.InterLocal, ProteinPropsMod.InterDistant]
    )
    push!(reg_sites, inter_inter)

    #prod_sites
    prod_sites = Array{ProteinProps, 1}()
    intra = rand_site(
        config,
        target=ProteinPropsMod.Intra
    )
    push!(prod_sites, intra)

    inter = rand_site(
        config,
        target=[ProteinPropsMod.InterLocal, ProteinPropsMod.InterDistant]
    )
    push!(prod_sites, inter)
    
    Gene(
        config,
        genome_index,
        reg_sites,
        prod_sites
    )
end

end

module GeneMod

using RunMod
using ProteinPropsMod
using MiscUtilsMod

import CustomEnumMod
import RandUtilsMod
import Base.show

export Gene

RegSites = CustomEnumMod.RegSites()
ProdSites = CustomEnumMod.ProdSites()

mutable struct Gene
    config::Config
    genome_index::Int64
    reg_sites::Array{ProteinProps, 1}
    prod_sites::Array{ProteinProps, 1}
end

function show(io::IO, gene::Gene, ilevel::Int64=0)
    iprintln(io, "Gene:", ilevel)
    iprintln(io, "genome_index: $(gene.genome_index)", ilevel + 1)
    
    iprintln(io, "reg_sites:", ilevel + 1)
    for (sym, val) in RegSites
        iprint(io, "$(string(sym)): ", ilevel + 2)
        ProteinPropsMod.show(io, gene.reg_sites[val], 0)
    end

    iprintln(io, "prod_sites:", ilevel + 1)
    for (sym, val) in instances(ProdSites)
        iprint(io, "$(string(sym)): ", ilevel + 2)
        ProteinPropsMod.show(io, gene.prod_sites[val], 0)
    end
end

function rand_site(
    config::Config;
    type::Union{ProteinPropsMod.ProteinType, Nothing}=nothing,
    target::Union{ProteinPropsMod.ProteinTarget, Nothing}=nothing,
    reg_action::Union{ProteinPropsMod.ProteinRegAction, Nothing}=nothing,
    app_action::Union{ProteinPropsMod.ProteinAppAction, Nothing}=nothing
)
    if type == nothing
        type = RandUtilsMod.rand_enum_val(config, ProteinPropsMod.ProteinType)
    end
    if target == nothing
        target = RandUtilsMod.rand_enum_val(config, ProteinPropsMod.ProteinTarget)
    end
    if reg_action == nothing
        reg_action = RandUtilsMod.rand_enum_val(config, ProteinPropsMod.ProteinRegAction)
    end
    if app_action == nothing
        app_action = RandUtilsMod.rand_enum_val(config, ProteinPropsMod.ProteinAppAction)
    end

    ProteinProps(type, target, reg_action, app_action)
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
        target=ProteinPropsMod.Inter
    )
    push!(reg_sites, inter_intra)

    inter_inter = rand_site(
        config,
        type=ProteinPropsMod.Reg,
        target=ProteinPropsMod.Inter
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
        target=ProteinPropsMod.Inter
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

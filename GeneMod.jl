module GeneMod

using RunMod
using ProteinMod

import RandUtilsMod
import Random

export Gene,
    rand_init

#accepts, regulates
@enum RegSites::Int64 IntraIntra=1 IntraInter=2 InterIntra=3 InterInter=4
#produces
@enum ProdSites::Int64 Intra=1 Inter=2

mutable struct Gene
    run::Run
    genome_index::Int64
    reg_sites::Array{ProteinProps, 1}
    prod_sites::Array{ProteinProps, 1}
end

function rand_site(
    run::Run;
    type::Union{ProteinMod.ProteinType, Nothing}=nothing,
    target::Union{ProteinMod.ProteinTarget, Nothing}=nothing,
    reg_action::Union{ProteinMod.ProteinRegAction, Nothing}=nothing
    app_action::Union{ProteinMod.ProteinAppAction, Nothing}=nothing
)
    if type == nothing
        type = RandUtilsMod.rand_enum_val(run, ProteinMod.ProteinType)
    end
    if target == nothing
        target = RandUtilsMod.rand_enum_val(run, ProteinMod.ProteinTarget)
    end
    if reg_action == nothing
        reg_action = RandUtilsMod.rand_enum_val(run, ProteinMod.ProteinRegAction)
    end
    if app_action == nothing
        app_action = RandUtilsMod.rand_enum_val(run, ProteinMod.ProteinAppAction)
    end

    ProteinProps(type, target, reg_action, app_action)
end

function rand_init(run::Run, genome_index::Int64)
    #reg_sites
    reg_sites = Array{ProteinProps, 1}()
    intra_intra = rand_site(
        run,
        type=ProteinMod.Reg,
        target=ProteinMod.Intra
    )
    push!(reg_sites, intra_intra)

    intra_inter = rand_site(
        run,
        type=ProteinMod.Reg,
        target=ProteinMod.Intra
    )
    push!(reg_sites, intra_inter)

    inter_intra = rand_site(
        run,
        type=ProteinMod.Reg,
        target=ProteinMod.Inter
    )
    push!(reg_sites, inter_intra)

    inter_inter = rand_site(
        run,
        type=ProteinMod.Reg,
        target=ProteinMod.Inter
    )
    push!(reg_sites, inter_inter)

    #prod_sites
    prod_sites = Array{ProteinProps, 1}()
    intra = rand_site(
        run,
        target=ProteinMod.Intra
    )
    push!(prod_sites, intra)

    inter = rand_site(
        run,
        target=ProteinMod.Inter
    )
    push!(prod_sites, inter)
    
    Gene(
        run,
        genome_index,
        reg_sites,
        prod_sites
    )
end

end

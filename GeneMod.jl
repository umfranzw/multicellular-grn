module GeneMod

using RunMod
using ProteinPropsMod

import RandUtilsMod

export Gene

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
    type::Union{ProteinPropsMod.ProteinType, Nothing}=nothing,
    target::Union{ProteinPropsMod.ProteinTarget, Nothing}=nothing,
    reg_action::Union{ProteinPropsMod.ProteinRegAction, Nothing}=nothing,
    app_action::Union{ProteinPropsMod.ProteinAppAction, Nothing}=nothing
)
    if type == nothing
        type = RandUtilsMod.rand_enum_val(run, ProteinPropsMod.ProteinType)
    end
    if target == nothing
        target = RandUtilsMod.rand_enum_val(run, ProteinPropsMod.ProteinTarget)
    end
    if reg_action == nothing
        reg_action = RandUtilsMod.rand_enum_val(run, ProteinPropsMod.ProteinRegAction)
    end
    if app_action == nothing
        app_action = RandUtilsMod.rand_enum_val(run, ProteinPropsMod.ProteinAppAction)
    end

    ProteinProps(type, target, reg_action, app_action)
end

function rand_init(run::Run, genome_index::Int64)
    #reg_sites
    reg_sites = Array{ProteinProps, 1}()
    intra_intra = rand_site(
        run,
        type=ProteinPropsMod.Reg,
        target=ProteinPropsMod.Intra
    )
    push!(reg_sites, intra_intra)

    intra_inter = rand_site(
        run,
        type=ProteinPropsMod.Reg,
        target=ProteinPropsMod.Intra
    )
    push!(reg_sites, intra_inter)

    inter_intra = rand_site(
        run,
        type=ProteinPropsMod.Reg,
        target=ProteinPropsMod.Inter
    )
    push!(reg_sites, inter_intra)

    inter_inter = rand_site(
        run,
        type=ProteinPropsMod.Reg,
        target=ProteinPropsMod.Inter
    )
    push!(reg_sites, inter_inter)

    #prod_sites
    prod_sites = Array{ProteinProps, 1}()
    intra = rand_site(
        run,
        target=ProteinPropsMod.Intra
    )
    push!(prod_sites, intra)

    inter = rand_site(
        run,
        target=ProteinPropsMod.Inter
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

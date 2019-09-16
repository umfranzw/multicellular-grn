module GeneMod

using RunMod
using BindSiteMod
using ProteinMod

import RandUtilsMod

export Gene,
    ProdLogic,
    rand_init

@enum ProdLogic::Bool ProdLogicAnd=false ProdLogicOr=true

mutable struct Gene
    run::Run
    genome_index::Int64
    initial_output_rate::Float64
    output_rate::Float64
    bind_threshold::Float64
    bind_sites::Array{BindSite, 1}
    prod_sites::Array{BindSite, 1}
    prod_logic::ProdLogic
    active_products::Array{Union{Protein, Nothing}, 1}
end

function rand_init(run::Run, genome_index::Int64)
    initial_output_rate = RandUtilsMod.rand_float(run)
    Gene(
        run,
        genome_index,
        initial_output_rate,
        initial_output_rate,
        RandUtilsMod.rand_float(run),
        map(x -> BindSiteMod.rand_init(run, ProteinMod.num_bits), 1:run.num_bind_sites),
        map(x -> BindSiteMod.rand_init(run, ProteinMod.num_bits), 1:run.num_bind_sites),
        RandUtilsMod.rand_enum_val(run, ProdLogic),
        []
    )
end

end

module GeneMod

using RunMod
using ProteinMod

import RandUtilsMod

export Gene,
    ProdLogic,
    rand_init

@enum ProdLogic::Bool ProdLogicAnd=false ProdLogicOr=true

mutable struct Gene
    run::Run
    genome_index::Int64
    bind_sites::Array{BitArray{1}, 1}
    prod_sites::Array{BitArray{1}, 1}
    prod_logic::ProdLogic
end

function rand_init(run::Run, genome_index::Int64)
    bind_sites = []
    prod_sites = []
    for i in 1:run.num_bind_sites
        push!(bind_sites, RandUtilsMod.rand_bits(run, ProteinMod.num_bits))
        push!(prod_sites, RandUtilsMod.rand_bits(run, ProteinMod.num_bits))
    end
    
    Gene(
        run,
        genome_index,
        bind_sites,
        prod_sites,
        RandUtilsMod.rand_enum_val(run, ProdLogic)
    )
end

end

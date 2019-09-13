module GeneMod

using RunMod
using BindSite
using ProteinMod

export Gene,
    ProdLogic

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

end

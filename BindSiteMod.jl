module BindSiteMod

using RunMod
using ProteinMod
import RandUtilsMod

export BindSite,
    rand_init

mutable struct BindSite
    seq::BitArray{1}
    bound_protein::Union{Protein, Nothing}
end

function rand_init(run::Run, n::Int64)
    BindSite(RandUtilsMod.rand_bits(run, n), nothing)
end

end

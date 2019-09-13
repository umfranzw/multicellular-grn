module BindSiteMod

using ProteinMod

export BindSite

mutable struct BindSite
    seq::BitArray{1}
    bound_protein::Union{Protein, Nothing}
end

end

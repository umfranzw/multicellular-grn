module CellMod

using GenomeMod

mutable struct Cell
    proteins::Dict{BitArray, Protein}
    genome::Genome
end

end

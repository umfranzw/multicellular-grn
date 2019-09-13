module CellMod

mutable struct Cell
    proteins::Dict{BitArray, Protein}
end

end

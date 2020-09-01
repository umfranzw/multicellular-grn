module AgeStepMod

using IndividualMod
using CellTreeMod

function run_age(indiv::Individual)
    CellTreeMod.traverse(cell -> cell.age += 1, indiv.cell_tree)
end

end

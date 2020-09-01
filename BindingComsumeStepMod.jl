module BindingComsumeStepMod

using IndividualMod
using CellTreeMod
using GeneStateMod
using CellMod

function run_binding_consum(indiv::Individual)
    CellTreeMod.traverse(run_binding_consum_for_cell, indiv.cell_tree)
end

function run_binding_consum_for_cell(cell::Cell)
    foreach(GeneStateMod.run_binding_consum, cell.gene_states)
end

end

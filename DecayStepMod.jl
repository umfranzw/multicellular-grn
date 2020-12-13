module DecayStepMod

using IndividualMod
using CellTreeMod
using CellMod
using ProteinStoreMod
using ProteinPropsMod
using GeneMod

function run_decay(indiv::Individual)
    CellTreeMod.traverse(run_decay_for_cell, indiv.cell_tree)
end

function run_decay_for_cell(cell::Cell)
    proteins = ProteinStoreMod.get_all(cell.proteins)
    for protein in proteins
        #decrease the concentration
        protein.concs = max.(protein.concs .- cell.config.run.decay_rate, zeros(length(protein.concs)))

        #remove any proteins that have decayed below the allowable threshold
        if all(c -> c < cell.config.run.protein_deletion_threshold, protein.concs)
            #clear any bindings that this protein has
            type = ProteinPropsMod.get_fcn(protein.props) == ProteinPropsMod.Inhibit ? ProdSite : BindSite
            for gs in cell.gene_states
                for i in 1:length(gs.bindings[type])
                    if gs.bindings[type][i] == protein
                        gs.bindings[type][i] = nothing
                    end
                end
            end

            #remove it from the cell
            ProteinStoreMod.remove(cell.proteins, protein)
        end
    end
end

end

module RegSimMod

using RunMod
using IndividualMod
using FitnessMod
using Printf
using CellTreeMod
import TrackerMod

function reg_sim(run::Run, pop::Array{Individual, 1}, ea_step::Int64)
    pop_trees = Array{Array{CellTree, 1}, 1}()
    for pop_index in 1:length(pop)
        #@info @sprintf("Individual %d\n", pop_index)
        indiv = pop[pop_index]
        indiv_trees = Array{CellTree, 1}()
        
        for reg_step in 1:run.reg_steps
            #@info @sprintf("Reg step %d\n", reg_step)

            IndividualMod.run_bind(indiv)
            IndividualMod.run_produce(indiv)
            IndividualMod.run_diffuse(indiv)
            IndividualMod.run_neighbour_comm(indiv)
            IndividualMod.run_protein_app(indiv)

            #update gene-protein interaction graph (used by mutation op)
            IndividualMod.update_chains(indiv, reg_step)
            #save before decay so we can reconstruct the concs and bindings properly later
            TrackerMod.save_reg_state(indiv.cell_tree, ea_step, reg_step, pop_index)
            IndividualMod.run_decay(indiv)
            IndividualMod.run_age(indiv)
        end

        IndividualMod.run_fix_syms(indiv)
        FitnessMod.eval(indiv, ea_step)
    end
end

end

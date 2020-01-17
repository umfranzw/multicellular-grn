module RegSimMod

using RunMod
using IndividualMod
using FitnessMod
using Printf
using CellTreeMod
import TrackerMod

reg_ops = (IndividualMod.run_bind,
           IndividualMod.run_produce,
           IndividualMod.run_diffuse,
           IndividualMod.run_protein_app,
           IndividualMod.run_decay)

function reg_sim(run::Run, pop::Array{Individual, 1}, ea_step::Int64)
    pop_trees = Array{Array{CellTree, 1}, 1}()
    for pop_index in 1:length(pop)
        #@info @sprintf("Individual %d\n", pop_index)
        indiv = pop[pop_index]
        indiv_trees = Array{CellTree, 1}()
        
        for reg_step in 1:run.reg_steps
            #@info @sprintf("Reg step %d\n", reg_step)

            for op in reg_ops
                op(indiv)
            end
            TrackerMod.save_reg_state(indiv.cell_tree, ea_step, reg_step, pop_index)
        end

        FitnessMod.eval(indiv)
    end
end

end

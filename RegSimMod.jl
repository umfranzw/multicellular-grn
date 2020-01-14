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

function reg_sim(run::Run, pop::Array{Individual, 1})
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
            push!(indiv_trees, deepcopy(indiv.cell_tree)) #!!!!
        end

        push!(pop_trees, indiv_trees)
        FitnessMod.eval(indiv)
    end

    TrackerMod.save_reg_state("after_reg_step", pop_trees)
end

end

module RegSimMod

using RunMod
using IndividualMod
using FitnessMod

import BindStepMod
import ProduceStepMod
import BindingComsumeStepMod
import DiffuseStepMod
import NeighbourCommStepMod
import ProteinAppStepMod
import DecayStepMod
import AgeStepMod
import TrackerMod

function reg_sim(run::Run, pop::Array{Individual, 1}, ea_step::Int64)
    if run.multithreaded
        Threads.@threads for pop_index in 1:length(pop)
            step(run, pop, pop_index, ea_step)
        end
    else
        for pop_index in 1:length(pop)
            step(run, pop, pop_index, ea_step)
        end
    end
end

#TODO:
#no need to save state at reg step 0 (but still need to save initial fitnesses on "ea_step 0")

function step(run::Run, pop::Array{Individual, 1}, pop_index::Int64, ea_step::Int64)
    #@info @sprintf("Individual %d\n", pop_index)
    indiv = pop[pop_index]
    #save initial state under index 0
    #TrackerMod.save_reg_state(indiv, ea_step, pop_index, 0, TrackerMod.IndivStateAfterBind)

    reg_step = 1
    while reg_step <= run.reg_steps && IndividualMod.has_proteins(indiv)
        #@info @sprintf("Reg step %d\n", reg_step)

        BindStepMod.run_bind(indiv)
        TrackerMod.save_reg_state(indiv, ea_step, pop_index, reg_step, TrackerMod.AfterBind)

        ProduceStepMod.run_produce(indiv)
        BindingComsumeStepMod.run_binding_consum(indiv)
        TrackerMod.save_reg_state(indiv, ea_step, pop_index, reg_step, TrackerMod.AfterProd)
        
        DiffuseStepMod.run_diffuse(indiv)
        NeighbourCommStepMod.run_neighbour_comm(indiv)
        ProteinAppStepMod.run_protein_app(indiv)
        DecayStepMod.run_decay(indiv)
        AgeStepMod.run_age(indiv)

        #FitnessMod.eval_intermediate(indiv, reg_step)

        reg_step += 1
    end

    IndividualMod.run_fix_syms(indiv)
    FitnessMod.eval_final(indiv, ea_step)
    #save final state under index run.reg_step + 1
    TrackerMod.save_reg_state(indiv, ea_step, pop_index, run.reg_steps + 1, TrackerMod.AfterBind)
end

end

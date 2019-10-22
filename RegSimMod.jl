module RegSimMod

using RunMod
using IndividualMod
using FitnessMod
using Printf

function reg_sim(run::Run, pop::Array{Individual, 1})
    for pop_index in 1:length(pop)
        #@info @sprintf("Individual %d\n", pop_index)
        indiv = pop[pop_index]
        
        for reg_step in 1:run.reg_steps
            #@info @sprintf("Reg step %d\n", reg_step)
            
            IndividualMod.run_bind(indiv)
            IndividualMod.run_produce(indiv)
            IndividualMod.run_diffuse(indiv)
            IndividualMod.run_protein_app(indiv)
            IndividualMod.run_decay(indiv)
        end

        FitnessMod.eval(indiv)
    end
end

end

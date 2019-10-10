module RegSimMod

using RunMod
using IndividualMod
using FitnessMod

function reg_sim(run::Run, pop::Array{Individual, 1})
    for indiv in pop
        for reg_step in 1:run.reg_steps
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

module EvAlgMod

using RunMod
using IndividualMod
using MutateMod
using RegSimMod

function ev_alg(run::Run)
    pop = map(i -> IndividualMod.rand_init(run), 1:run.pop_size)

    for ea_step in 1:run.ea_steps
        @info "EA step" ea_step
        
        #run the genetic operator
        MutateMod.mutate(run, pop)
        
        #the reg sim will update the fitnesses
        RegSimMod.reg_sim(run, pop)

        #reset the individuals before the next iteration
        map(IndividualMod.reset, pop)
    end
end

end

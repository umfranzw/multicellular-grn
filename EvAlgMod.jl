module EvAlgMod

using RunMod
using IndividualMod
using MutateMod
using RegSimMod
using TrackerMod
using Printf

function ev_alg(run::Run)
    TrackerMod.create_tracker(run)
    pop = map(i -> IndividualMod.rand_init(run), 1:run.pop_size)

    for ea_step in 1:run.ea_steps
        @info @sprintf("EA step: %d", ea_step)
        TrackerMod.save_pop_state(ea_step, pop)
        
        #run the genetic operator
        MutateMod.mutate(pop)
        
        #the reg sim will update the fitnesses
        RegSimMod.reg_sim(run, pop)

        TrackerMod.update_bests(ea_step, pop)

        #reset the individuals before the next iteration
        map(IndividualMod.reset, pop)
    end

    TrackerMod.destroy_tracker()
end

end

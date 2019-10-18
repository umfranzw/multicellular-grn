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
        
        #run the genetic operator
        MutateMod.mutate(run, pop)
        
        #the reg sim will update the fitnesses
        RegSimMod.reg_sim(run, pop)
        TrackerMod.track_state(ea_step, pop)

        #reset the individuals before the next iteration
        map(IndividualMod.reset, pop)
    end

    TrackerMod.destroy_tracker()
end

end

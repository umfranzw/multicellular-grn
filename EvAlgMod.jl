module EvAlgMod

using RunMod
using IndividualMod
using MutateMod
using RegSimMod
using TrackerMod
using Printf
import Random

gen_ops = (MutateMod.mutate)

function ev_alg(run::Run)
    TrackerMod.create_tracker(run)
    pop = create_pop(run)
    TrackerMod.save_state("initial_pop", pop)
    RegSimMod.reg_sim(run, pop)
    TrackerMod.update_bests(pop)

    ea_step = 1
    while !terminate(run) && ea_step < run.ea_steps
        @info @sprintf("EA step: %d", ea_step)
        
        #run the genetic operators
        for op in gen_ops
            op(pop)
        end
        TrackerMod.save_state_on_step(ea_step, "after_gen_ops", pop)
        
        #the reg sim will update the fitnesses
        RegSimMod.reg_sim(run, pop)

        TrackerMod.update_bests(ea_step, pop)

        #reset the individuals before the next iteration
        map(IndividualMod.reset, pop)
    end

    TrackerMod.destroy_tracker()
end

function terminate(run::Run)
    TrackerMod.run_best != nothing && TrackerMod.run_best.fitness <= run.fitness_term_threshold
end

function create_pop(run::Run)
    pop = Array{Individual, 1}()
    if run.fix_rng_seed
        seed_base = run.rng_seed
    else
        dev = Random.RandomDevice()
        seed_base = UInt64(Random.rand(dev) * 0xffffffffffffffff)
    end
    
    for i in 1:run.pop_size
        push!(pop, IndividualMod.rand_init(run, seed_base + i))
    end

    pop
end

end

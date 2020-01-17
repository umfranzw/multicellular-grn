module EvAlgMod

using RunMod
using IndividualMod
using MutateMod
using RegSimMod
using TrackerMod
using Printf
import Random

gen_ops = (MutateMod.mutate,)

function ev_alg(run::Run)
    data_filename = join((RunMod.DATA_PATH, run.data_output_file), "/")
    TrackerMod.create_tracker(run, data_filename)
    pop = create_pop(run)
    #TrackerMod.save_ea_state("initial_pop", pop)
    RegSimMod.reg_sim(run, pop, -1)
    TrackerMod.update_bests(pop)

    ea_step = 1
    while !terminate(run) && ea_step <= run.ea_steps
        @info @sprintf("EA step: %d", ea_step)
        
        #run the genetic operators
        for op in gen_ops
            op(pop)
        end
        TrackerMod.save_ea_state(pop, ea_step)
        
        #the reg sim will update the fitnesses
        RegSimMod.reg_sim(run, pop, ea_step)

        TrackerMod.update_bests(pop)

        #reset the individuals before the next iteration
        map(IndividualMod.reset, pop)

        ea_step += 1
    end

    TrackerMod.destroy_tracker()
end

function terminate(run::Run)
    tracker = TrackerMod.get_tracker()
    
    tracker.run_best != nothing && tracker.run_best.fitness <= run.fitness_term_threshold
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

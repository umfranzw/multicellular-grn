module EvAlgMod

using RunMod
using IndividualMod
using MutateMod
using SelectionMod
using RegSimMod
using TrackerMod
using BestInfoMod
using Printf
import Random

function ev_alg(run::Run)
    data_filename = join((RunMod.DATA_PATH, run.data_output_file), "/")
    pop = create_pop(run)
    selector = Selector(run)
    TrackerMod.create_tracker(run, data_filename)
    #TrackerMod.save_ea_state(pop, 0, true)
    RegSimMod.reg_sim(run, pop, 0)
    TrackerMod.update_fitnesses(pop, 0)
    foreach(IndividualMod.reset_cell_tree, pop)

    ea_step = 1
    while !terminate(run) && ea_step <= run.ea_steps
        @info @sprintf("EA step: %d", ea_step)
        
        #run the genetic operators
        pop = SelectionMod.select(selector, pop)
        MutateMod.mutate(run, pop, ea_step)

        #TrackerMod.save_ea_state(pop, ea_step)
        foreach(IndividualMod.reset_reg_sim_info, pop)
        
        #the reg sim will update the fitnesses
        RegSimMod.reg_sim(run, pop, ea_step)

        TrackerMod.update_fitnesses(pop, ea_step)
        #print_trees(pop)

        #reset the individuals before the next iteration
        foreach(IndividualMod.reset_cell_tree, pop)

        ea_step += 1
    end

    TrackerMod.save_run_best()
    TrackerMod.save_fitnesses()
    TrackerMod.destroy_tracker()
end

function print_trees(pop::Array{Individual, 1})
    foreach(i -> println(i.cell_tree), pop)
    print("\n")
end

function print_fitnesses(pop::Array{Individual, 1})
    foreach(i -> println(@sprintf("%0.2f", i.fitness)), pop)
    print("\n")
end

function print_pop(pop::Array{Individual, 1})
    for indiv in pop
        println(indiv)
    end

    print("\n")
end

function terminate(run::Run)
    tracker = TrackerMod.get_tracker()
    
    BestInfoMod.is_set(tracker.run_best) && tracker.run_best.indiv.fitness <= run.fitness_term_threshold
end

function create_pop(run::Run)
    pop = Array{Individual, 1}()
    for i in 1:run.pop_size
        push!(pop, IndividualMod.rand_init(run, UInt64(i)))
    end

    pop
end

end

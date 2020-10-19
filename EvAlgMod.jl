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
    @info "Setting up"
    data_filename = join((RunMod.DATA_PATH, run.data_output_file), "/")
    pop = create_pop(run)
    selector = Selector(run)
    TrackerMod.create_tracker(run, data_filename)
    #TrackerMod.save_ea_state(pop, 0, true)

    @info "Beginning algorithm"

    @time begin
        @info "EA step: 0"
        RegSimMod.reg_sim(run, pop, 0)
        TrackerMod.update_fitnesses(pop, 0)
        foreach(IndividualMod.reset_cell_tree, pop)

        ea_step = 1
        while !terminate(run) && ea_step <= run.ea_steps
            @info @sprintf("EA step: %d", ea_step)
            
            #run the genetic operators
            pop = SelectionMod.select(selector, pop)
            prev_pop = deepcopy(pop)
            MutateMod.mutate(run, pop, ea_step)

            #TrackerMod.save_ea_state(pop, ea_step)
            foreach(IndividualMod.reset_reg_sim_info, pop)
            
            #the reg sim will update the fitnesses in the indiv objects
            RegSimMod.reg_sim(run, pop, ea_step)

            #don't allow individuals whose fitness got worse because of mutation into the next generation (keep the indiv from prev_pop in this case)
            enforce_fitness_front(run, prev_pop, pop)
            
            #update fitnesses in the fitnesses array, as well as the best fitnesses
            TrackerMod.update_fitnesses(pop, ea_step)
            #print_trees(pop)

            #reset the individuals before the next iteration
            foreach(IndividualMod.reset_cell_tree, pop)

            ea_step += 1
        end
    end

    @info "Cleaning up"
    TrackerMod.save_run_best()
    TrackerMod.save_fitnesses()
    TrackerMod.destroy_tracker()
    @info "Done"
end

function enforce_fitness_front(run::Run, prev_pop::Array{Individual, 1}, pop::Array{Individual, 1})
    if run.multithreaded
        Threads.@threads for i in 1:length(prev_pop)
            if prev_pop[i].fitness < pop[i].fitness
                pop[i] = prev_pop[i]
            end
        end

    else
        for i in 1:length(prev_pop)
            if prev_pop[i].fitness < pop[i].fitness
                pop[i] = prev_pop[i]
            end
        end
    end
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
    pop = Array{Individual, 1}(undef, run.pop_size)
    if run.multithreaded
        Threads.@threads for i in 1:run.pop_size
            pop[i] = IndividualMod.rand_init(run, UInt64(i))
        end
        
    else
        for i in 1:run.pop_size
            pop[i] = IndividualMod.rand_init(run, UInt64(i))
        end
    end

    pop
end

end

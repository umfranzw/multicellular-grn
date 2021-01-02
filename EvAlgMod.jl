module EvAlgMod

using RunMod
using IndividualMod
using CrossoverMod
using MutateMod
using SelectionMod
using RegSimMod
using TrackerMod
using BestInfoMod
using GrowthMod
using Printf
import MiscUtilsMod
import Random

function ev_alg(run::Run)
    @info "Setting up"
    data_filename = join((RunMod.DATA_PATH, run.data_output_file), "/")
    pop = create_pop(run)
    TrackerMod.create_tracker(run, data_filename)
    #TrackerMod.save_ea_state(pop, 0, true)

    @info "Beginning algorithm"

    val, el_time, total_bytes, gc_time, counters = @timed begin
        step_output_buf = IOBuffer()
        write_ea_step_title(step_output_buf, 0)
        RegSimMod.reg_sim(run, pop, 0)
        TrackerMod.update_fitnesses(pop, 0, step_output_buf)
        foreach(IndividualMod.reset_cell_tree, pop)
        @info String(take!(step_output_buf))
        pop = SelectionMod.select(run, pop)

        ea_step = 1
        while !terminate(run) && ea_step <= run.ea_steps
            write_ea_step_title(step_output_buf, ea_step)

            #save old pop
            #prev_pop = deepcopy(pop)
            
            #do growth
            GrowthMod.grow(run, pop)
            
            #do mutation
            MutateMod.mutate(run, pop, ea_step)

            #do crossover
            CrossoverMod.crossover(run, pop, ea_step)

            #TrackerMod.save_ea_state(pop, ea_step)
            foreach(IndividualMod.reset_reg_sim_info, pop)
            
            #the reg sim will update the fitnesses in the indiv objects
            RegSimMod.reg_sim(run, pop, ea_step)

            #don't allow individuals whose fitness got worse because of genetic operations into the next generation (keep the indiv from prev_pop in this case)
            #enforce_fitness_front(run, prev_pop, pop, ea_step)
            
            #update fitnesses in the fitnesses array, as well as the best fitnesses
            TrackerMod.update_fitnesses(pop, ea_step, step_output_buf)
            #print_trees(pop)

            #reset the individuals before the next iteration
            foreach(IndividualMod.reset_cell_tree, pop)

            #do selection for next generation
            pop = SelectionMod.select(run, pop)

            ea_step += 1
            
            @info String(take!(step_output_buf))
        end
    end

    @info "Logging final data"
    TrackerMod.save_run_best()
    TrackerMod.save_fitnesses()
    data_file_size = TrackerMod.destroy_tracker() / 2^20 #in MiB
    total_alloc = total_bytes / 2^20 #in MiB
    num_allocs = Base.gc_alloc_count(counters)
    
    write(step_output_buf, "\n------\n")
    write(step_output_buf, "Stats:\n")
    write(step_output_buf, "------\n")
    write(step_output_buf, @sprintf("Elapsed time: %s (%s gc time)\n", MiscUtilsMod.get_time_str(el_time), MiscUtilsMod.get_time_str(gc_time)))
    write(step_output_buf, @sprintf("Memory: %0.2f MiB (%d allocations)\n", total_alloc, num_allocs))
    if run.log_level > RunMod.LogNone
        write(step_output_buf, @sprintf("Data file size: %0.2f MiB", data_file_size))
    end
    
    @info String(take!(step_output_buf))
    @info "Done"
end

function write_ea_step_title(output_buf::IOBuffer, ea_step::Int64)
    title = "EA step: $(ea_step)"
    divider = repeat('-', length(title))
    write(output_buf, "\n$(divider)\n")
    write(output_buf, "$(title)\n")
    write(output_buf, "$(divider)\n")
end

function enforce_fitness_front(run::Run, prev_pop::Array{Individual, 1}, pop::Array{Individual, 1}, ea_step::Int64)
    if run.multithreaded
        Threads.@threads for i in 1:length(prev_pop)
            if prev_pop[i].fitness < pop[i].fitness
                pop[i] = prev_pop[i]
                TrackerMod.save_rolledback_state(ea_step, i)
            end
        end

    else
        for i in 1:length(prev_pop)
            if prev_pop[i].fitness < pop[i].fitness
                pop[i] = prev_pop[i]
                TrackerMod.save_rolledback_state(ea_step, i)
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

module SelectionMod

using IndividualMod
using RunMod
using RandUtilsMod
import Random

function select(run::Run, pop::Array{Individual, 1})
    new_pop = Array{Individual, 1}(undef, run.pop_size)

    #try to maintain order for the sake of making vis analysis easier
    #Guarentee: if an indiv was not completely weeded out of the pop, it will be in its original position.
    #Identical copies of the indiv may also occupy the positions of other indivs that have been weeded out.
    old_order = Array{UInt64, 1}(undef, run.pop_size)
    new_order = Dict{UInt64, Array{Int64, 1}}()

    if run.multithreaded
        Threads.@threads for i in 1:run.pop_size
            old_order[i] = pop[i].id
            
            tourn = Random.rand(pop[i].config.rng, pop, run.tourn_size)
            winner = foldl((cur_best, next) -> next.fitness < cur_best.fitness ? next : cur_best, tourn)
            #insert a *copy*, since we may pick the same winner more than once
            new_pop[i] = deepcopy(winner)
        end

        for i in 1:run.pop_size
            if new_pop[i].id in keys(new_order)
                push!(new_order[new_pop[i].id], i)
            else
                new_order[new_pop[i].id] = Array{Int64, 1}([i])
            end
        end
    else
        for i in 1:run.pop_size
            old_order[i] = pop[i].id
            
            tourn = Random.rand(pop[i].config.rng, pop, run.tourn_size)
            winner = foldl((cur_best, next) -> next.fitness < cur_best.fitness ? next : cur_best, tourn)
            #insert a *copy*, since we may pick the same winner more than once
            new_pop[i] = deepcopy(winner)

            if new_pop[i].id in keys(new_order)
                push!(new_order[new_pop[i].id], i)
            else
                new_order[new_pop[i].id] = Array{Int64, 1}([i])
            end
        end
    end

    for i in 1:run.pop_size
        old_id = old_order[i]
        new_id = new_pop[i].id
        if old_id != new_id && old_id in keys(new_order)
            #prefer indivs with matching last_mod attributes if possible
            index = findfirst(new_index -> new_pop[new_index].last_mod == pop[i].last_mod, new_order[old_id])
            #otherwise just take the back one
            if index == nothing
                swap_index = pop!(new_order[old_id])
            else
                swap_index = new_order[old_id][index]
                deleteat!(new_order[old_id], index)
            end
            
            if length(new_order[old_id]) == 0
                delete!(new_order, old_id)
            end
            new_pop[i], new_pop[swap_index] = new_pop[swap_index], new_pop[i]
        end
    end
    
    new_pop
end

end

module SelectionMod

using IndividualMod
using RunMod
using RandUtilsMod
import Random

function select(run::Run, pop::Array{Individual, 1})
    new_pop = Array{Individual, 1}(undef, run.pop_size)

    if run.multithreaded
        Threads.@threads for i in 1:run.pop_size
            tourn = Random.rand(pop[i].config.rng, pop, run.tourn_size)
            winner = foldl((cur_best, next) -> next.fitness < cur_best.fitness ? next : cur_best, tourn)
            #insert a *copy*, since we may pick the same winner more than once
            new_pop[i] = deepcopy(winner)
        end
    else
        for i in 1:run.pop_size
            tourn = Random.rand(pop[i].config.rng, pop, run.tourn_size)
            winner = foldl((cur_best, next) -> next.fitness < cur_best.fitness ? next : cur_best, tourn)
            #insert a *copy*, since we may pick the same winner more than once
            new_pop[i] = deepcopy(winner)
        end
    end
    
    new_pop
end

end

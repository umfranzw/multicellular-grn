module SelectionMod

using IndividualMod
using RunMod
using RandUtilsMod

import Random

export Selector

mutable struct Selector
    config::Config

    function Selector(run::Run)
        if run.fix_rng_seed
            seed_base = run.rng_seed
        else
            dev = Random.RandomDevice()
            seed_base = UInt64(Random.rand(dev) * 0xffffffffffffffff)
        end

        seed = seed_base + run.pop_size + 1 #note: individuals generate configs that use offsets 1:run.pop_size
        rng = Random.MersenneTwister(seed)
        config = Config(run, rng)

        new(config)
    end
end

function select(selector::Selector, pop::Array{Individual, 1})
    num_tourns = selector.config.run.pop_size ÷ selector.config.run.tourn_size
    #get the indices of the sorted population (from worst (highest) to best (lowest))
    sorted_indices = sortperm(pop, by=indiv -> indiv.fitness, rev=true)

    replacements = Array{Individual, 1}()
    if selector.config.run.use_tourn_replacement #random selection (individuals may be picked twice
        for i in 1:num_tourns
            winner = nothing
            for j in 1:selector.config.run.tourn_size
                indiv = Random.rand(config.rng, pop)
                if winner == nothing || indiv.fitness < winner.fitness
                    winner = indiv
                end
            end
            push!(replacements, winner)
        end
        
    else #individuals should never be picked more than once
        perm = Random.shuffle(selector.config.rng, 1:selector.config.run.pop_size)
        for i in 1:num_tourns
            winner = nothing
            for j in 1:selector.config.run.tourn_size
                perm_index = (i - 1) * selector.config.run.tourn_size + j
                pop_index = perm[perm_index]
                indiv = pop[pop_index]
                if winner == nothing || indiv.fitness < winner.fitness
                    winner = indiv
                end
            end
            push!(replacements, winner)
        end
    end

    for i in 1:num_tourns
        replace_index = sorted_indices[i]
        #println("Replacing $(pop[replace_index].fitness) with $(replacements[i].fitness)")
        replacement = replacements[i]
        pop[replace_index] = deepcopy(replacements[i])
    end
end

end

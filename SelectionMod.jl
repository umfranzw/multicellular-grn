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
    new_pop = Array{Individual, 1}()
    for i in 1:selector.config.run.pop_size
        tourn = Random.rand(selector.config.rng, pop, selector.config.run.tourn_size)
        winner = foldl((cur_best, next) -> next.fitness < cur_best.fitness ? next : cur_best, tourn)
        #insert a *copy*, since we may pick the same winner more than once
        push!(new_pop, deepcopy(winner))
    end
    
    new_pop
end

end

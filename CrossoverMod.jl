module CrossoverMod

#Assumes selection has been done beforehand (pop contains the newly selected individuals)

using RunMod
using IndividualMod
using GeneMod
import Random

config = nothing

function crossover(run::Run, pop::Array{Individual, 1})
    global config
    
    if config == nothing
        config = Config(run, UInt64(length(pop) + 1))
    end

    #TODO: how many individuals do we want to cross?
    #use crossover_proportion param (eg. 40% of pop): pop_size * 0.4 rounded to the nearest even number
    
    p1_index = Random.rand(config.rng, 1:length(pop))
    p2_index = Random.rand(config.rng, 1:length(pop))
    p1 = pop[p1_index]
    p2 = pop[ps_index]

    split_p1 = Random.rand(config.rng, 1:length(p1.genes))
    split_p2 = Random.rand(config.rng, 1:length(p2.genes))

    c1 = make_child(p1, p2, 1:split_p1, split_p2:length(p2.genes))
    c2 = make_child(p2, p1, 1:split_p2, split_p1:length(p1.genes))

    pop[p1_index] = c1
    pop[p2_index] = c2
end

function make_child(p1::Indiv, p2::Indiv, p1_range::UnitRange{Int64}, p2_range::Int64{In64})
    #assumption: child will replace p1, so it can take over it's config
    #assumption: parents will both be replaced, so we can steal their genes without copying them
    genes = Array{Gene, 1}()
    genome_index = 1
    for gene in p1.genes[p1_range]
        gene.genome_index = genome_index
        push!(genes, gene)
        genome_index += 1
    end

    for gene in p2.genes[p2_range]
        gene.genome_index = genome_index
        push!(genes, gene)
        genome_index += 1
    end
    
    Individual(p1.config, genes)
end

end

module CrossoverMod

#Assumes selection has been done beforehand (pop contains the newly selected individuals)

using RunMod
using IndividualMod
using GeneMod
import Random

config = nothing

function crossover(run::Run, pop::Array{Individual, 1}, ea_step::Int64)
    global config
    
    if config == nothing
        config = Config(run, UInt64(length(pop) + 1))
    end

    #use crossover_proportion param (eg. 40% of pop): pop_size * 0.4 rounded to the nearest even number
    num_cross = Int64(floor(run.pop_size * run.cross_prop))
    num_pairs = num_cross ÷ 2
    #note: duplication of fitter individuals has already occurred in the selection
    #So here, we only allow each individual to participate in crossover once
    #This also ensures that children are not crossed on iterations after the one they were generated on
    p1_choices = Random.shuffle(config.rng, 1:length(pop))
    p2_choices = Random.shuffle(config.rng, 1:length(pop))
    for i in 1:num_pairs
        p1_index = pop!(p1_choices)
        p2_index = pop!(p2_choices)
        p1 = pop[p1_index]
        p2 = pop[p2_index]

        split_p1 = Random.rand(config.rng, 1:length(p1.genes))
        split_p2 = Random.rand(config.rng, 1:length(p2.genes))

        pop_seed_offset = run.pop_size + 1
        ea_step_seed_offset = ea_step * num_cross
        c1_offset = (i - 1) * 2 + 1
        c1 = make_child(p1, p2, 1:split_p1 - 1, split_p2:length(p2.genes), pop_seed_offset + ea_step_seed_offset + c1_offset)
        #check("c1 after c1", c1)

        c2_offset = i * 2 #note: this is equivalent to (i - 1) * 2 + 2
        c2 = make_child(p2, p1, 1:split_p2 - 1, split_p1:length(p1.genes), pop_seed_offset + ea_step_seed_offset + c2_offset)
        #check("c1 after c2", c1)
        #check("c2 after c2", c2)

        pop[p1_index] = c1
        pop[p2_index] = c2
    end
end

#for debugging
# function check(msg::String, child::Individual)
#     gs_indices = map(g -> g.gene.genome_index, child.cell_tree.root.gene_states)
#     if gs_indices != collect(1:length(gs_indices))
#         println(msg)
#         println("Gene States out")
#         println(gs_indices)
#         println()
#     end

#     g_indices = map(g -> g.genome_index, child.genes)
#     if g_indices != collect(1:length(g_indices))
#         println(msg)
#         println("Genes out")
#         println(g_indices)
#         println()
#     end
# end

function make_child(p1::Individual, p2::Individual, p1_range::UnitRange{Int64}, p2_range::UnitRange{Int64}, seed_offset::Int64)
    #note: although both children are replaced in the current population, it's possible that one of them may be rolled back to
    #its parent at the end of this ea_step (see EvAlgMod.enforce_fitness_front()). So we can't just steal the genes and config
    #from the parents - we have to deepcopy them.
    genes = Array{Gene, 1}()
    genome_index = 1
    for gene in p1.genes[p1_range]
        child_gene = deepcopy(gene)
        child_gene.genome_index = genome_index
        # reference = map(objectid, genes)
        # i = findfirst(g -> objectid(g) in reference, genes)
        # if i != nothing
        #     println("Ahoy there!1")
        # end
        push!(genes, child_gene)
        genome_index += 1
    end

    for gene in p2.genes[p2_range]
        child_gene = deepcopy(gene)
        child_gene.genome_index = genome_index
        # reference = map(objectid, genes)
        # i = findfirst(g -> objectid(g) in reference, genes)
        # if i != nothing
        #     println("Ahoy there!1")
        # end
        push!(genes, child_gene)
        genome_index += 1
    end

    config = Config(p1.config.run, UInt64(p1.config.rng.seed[1]) + UInt64(seed_offset))
    Individual(config, genes)
end

end

module MutateMod

using RunMod
using IndividualMod
using GeneMod
using GeneStateMod
using ProteinMod
using ProteinPropsMod
using ProteinStoreMod
using RandUtilsMod
using SettingsMod
using CellTreeMod
using RegSimInfoMod

import Random
import RandUtilsMod

function mutate(run::Run, pop::Array{Individual, 1}, ea_step::Int64)
    if run.multithreaded
        Threads.@threads for indiv in pop
            mutate_indiv(indiv, ea_step)
        end
    else
        for indiv in pop
            mutate_indiv(indiv, ea_step)
        end
    end
end

function mutate_indiv(indiv::Individual, ea_step::Int64)
    mutated = point_mutate_indiv(indiv)

    if mutated
        indiv.last_mod = ea_step
    end

    mutated
end

function point_mutate_indiv(indiv::Individual)
    mutated = false
    for gene in indiv.genes
        mutated = mutated || point_mutate_gene(indiv, gene)
    end

    mutated
end

function point_mutate_gene(indiv::Individual, gene::Gene)
    mutated = false

    #mutate prod site args
    for prod_site in gene.prod_sites
        if RandUtilsMod.rand_decision(indiv.config, indiv.config.run.point_mut_prob)
            sign = Int8(RandUtilsMod.rand_decision(indiv.config, 0.5) * 2 - 1)
            prod_site.arg = sign * Random.rand(indiv.config.rng, Int8(0):Int8(indiv.config.run.max_protein_arg))
            mutated = true
        end
    end

    #mutate initial protein args
    for initial_protein in indiv.initial_cell_proteins
        if RandUtilsMod.rand_decision(indiv.config, indiv.config.run.point_mut_prob)
            sign = Int8(RandUtilsMod.rand_decision(indiv.config, 0.5) * 2 - 1)
            new_arg = sign * Random.rand(indiv.config.rng, Int8(0):Int8(indiv.config.run.max_protein_arg))
            active_protein = ProteinStoreMod.get(indiv.cell_tree.root.proteins, initial_protein.props)
            initial_protein.props.arg = new_arg
            active_protein.props.arg = new_arg
            
            mutated = true
        end
    end

    mutated
end

end

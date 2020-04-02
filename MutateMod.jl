module MutateMod

using RunMod
using IndividualMod
using GeneMod
using ProteinMod
using ProteinPropsMod
using RandUtilsMod

import Random

function mutate(pop::Array{Individual, 1}, ea_step::Int64)
    for indiv in pop
        mutate_indiv(indiv, ea_step)
    end
end

function mutate_indiv(indiv::Individual, ea_step::Int64)
    score_total = sum(indiv.gene_scores)
    copy_locations = Array{Int64, 1}()
    gene_index = 1
    while gene_index <= length(indiv.genes)
        normal_score = indiv.gene_scores[gene_index] / score_total
        if normal_score > indiv.config.run.gene_score_threshold
            mut_copy = dup_and_mutate_gene(indiv.genes[gene_index], ea_step)
            if mut_copy != nothing
                insert!(indiv.genes, gene_index + 1, mut_copy) #insert mutated copy after the src gene
                insert!(indiv.gene_scores, gene_index + 1, 0) #insert a new score for the new gene

                #insert a new conc into the initial proteins
                for protein in indiv.initial_cell_proteins
                    insert!(protein.concs, gene_index + 1, RandUtilsMod.rand_float(indiv.config))
                end
                
                gene_index += 1 #skip over the copy
            end
        else
            point_mutate_gene(indiv.genes[gene_index], ea_step)
        end
        
        gene_index += 1
    end
end

function dup_and_mutate_gene(gene::Gene, ea_step::Int64)
    copy = nothing
    if RandUtilsMod.rand_float(gene.config) < gene.config.run.mut_prob
        copy = deepcopy(gene)
        point_mutate_gene(copy, ea_step)
    end

    copy
end

function point_mutate_gene(gene::Gene, ea_step::Int64)
    #bind sites:
    for site in gene.bind_sites
        mutate_bind_site(gene.config, site, ea_step)
    end

    #prod sites:
    for site in gene.prod_sites
        mutate_props(gene.config, site, ea_step)
    end
end

function mutate_bind_site(config::Config, site::BindSite, ea_step::Int64)
    if RandUtilsMod.rand_float(config) < config.run.mut_prob
        #type can be anything but Application
        site.type = Random.rand(config.rng, [ProteinPropsMod.Internal, ProteinPropsMod.Neighbour, ProteinPropsMod.Diffusion])
    end
    if RandUtilsMod.rand_float(config) < config.run.mut_prob
        site.action = RandUtilsMod.rand_enum_val(config, ProteinPropsMod.ProteinAction)
    end
    if RandUtilsMod.rand_float(config) < config.run.mut_prob
        site.loc = RandUtilsMod.rand_enum_val(config, ProteinPropsMod.ProteinLoc)
    end

    #mutate site.threshold and site.consum_rate
    mutate_floats(config, site, ea_step)
end

function mutate_floats(config::Config, site::BindSite, ea_step::Int64)
    time_factor = 1.0 - ea_step / config.run.ea_steps
    
    for fieldname in (:threshold, :consum_rate)
        if RandUtilsMod.rand_float(config) < config.run.mut_prob
            range = time_factor
            delta = RandUtilsMod.rand_float(config) * range - range / 2 #value in [-range / 2, +range / 2]
            cur_val = getfield(site, fieldname)
            new_val = clamp(cur_val + delta, 0.0, 1.0)
            setfield!(site, fieldname, new_val)
        end
    end
end

#note: just because a mutation happens doesn't necessarily mean anything has changed - it's possible for random mutations to choose the same val that currently exists...
function mutate_props(
    config::Config,
    props::ProteinProps,
    ea_step::Int64;
    type::Union{Array{ProteinPropsMod.ProteinType, 1}, Nothing}=nothing,
    loc::Union{Array{ProteinPropsMod.ProteinLoc, 1}, Nothing}=nothing,
    action::Union{Array{ProteinPropsMod.ProteinAction, 1}, Nothing}=nothing,
    arg::Union{Array{Float64, 1}, Nothing}=nothing
)
    mutated = false
    i = 1
    enum_info = (
        (ProteinPropsMod.ProteinType, type, :type),
        (ProteinPropsMod.ProteinLoc, loc, :loc),
        (ProteinPropsMod.ProteinAction, action, :action)
    )
    while !mutated && i <= length(enum_info)
        if RandUtilsMod.rand_float(config) < config.run.mut_prob
            enum, options, fieldname = enum_info[i]
            if options == nothing
                new_val = RandUtilsMod.rand_enum_val(config, enum)
            else
                new_val = Random.rand(config.rng, options)
            end
            setfield!(props, fieldname, new_val)
            mutated = true
        end
        i += 1
    end

    #mutate props.arg
    if !mutated && RandUtilsMod.rand_float(config) < config.run.mut_prob
        if arg == nothing
            max_int8 = 2^7 - 1
            time_factor = 1.0 - ea_step / config.run.ea_steps
            range = Int64(floor(time_factor * max_int8))
            delta = RandUtilsMod.rand_int(config, 0, range) - range ÷ 2
            props.arg = Int8(clamp(props.arg + delta, -max_int8, max_int8))
        else
            props.arg = Random.rand(config.rng, arg)
        end
        mutated = true
    end

    mutated
end

end

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

import Random

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
    # score_total = sum(indiv.gene_scores)
    # copy_locations = Array{Int64, 1}()
    gene_index = 1
    while gene_index <= length(indiv.genes)
        # if ea_step < indiv.config.run.gene_dup_gen_limit && indiv.gene_scores[gene_index] > indiv.config.run.gene_score_threshold && length(indiv.genes) < indiv.config.run.max_genes
        #     mut_copy = dup_and_mutate_gene(indiv.genes[gene_index], ea_step)
        #     if mut_copy != nothing
        #         insert!(indiv.genes, gene_index + 1, mut_copy) #insert mutated copy after the src gene
        #         insert!(indiv.cell_tree.root.gene_states, gene_index + 1, GeneState(indiv.config.run, mut_copy)) #insert a new GeneState into the root cell
        #         insert!(indiv.gene_scores, gene_index + 1, 0) #insert a new score for the new gene
        #         IndividualMod.extend_sensors(indiv, gene_index + 1)

        #         new_num_concs = length(indiv.genes)
        #         #insert a new conc into the initial proteins, and into the corresponding copy that is in the root cell
        #         for init_protein in indiv.initial_cell_proteins
        #             active_protein = ProteinStoreMod.get(indiv.cell_tree.root.proteins, init_protein.props)
        #             conc = RandUtilsMod.rand_float(indiv.config)
        #             insert!(init_protein.concs, gene_index + 1, conc)
                    
        #             #this is needed since there can be duplicate proteins with identical props in cell.initial_cell_proteins
        #             if length(active_protein.concs) < new_num_concs
        #                 insert!(active_protein.concs, gene_index + 1, conc)
        #             end
        #         end
                
        #         gene_index += 1 #skip over the copy
        #     end
        # else
            point_mutate_gene(indiv.genes[gene_index], ea_step)
        # end
        
        gene_index += 1
    end

    #swap locations
    mutate_location(indiv.config, indiv.genes)

    #update genome_indices
    # foreach(i -> indiv.genes[i].genome_index = i, 1:length(indiv.genes))
end

function mutate_location(config::Config, genes::Array{Gene, 1})
    if RandUtilsMod.rand_float(gene.config) < gene.config.run.mut_prob
        src = RandUtilsMod.rand_int(indiv.config, 1, length(genes))
        delta = Random.rand(indiv.config, (1, -1))
        dest = (src + delta) < 1 ? (src + delta) + length(genes) : (src + delta) % length(genes)
        genes[src], genes[dest] = genes[dest], genes[src]
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
        mutate_prod_site(gene.config, site, ea_step)
    end
end

function mutate_prod_site(config::Config, site::ProdSite, ea_step::Int64)
    if RandUtilsMod.rand_float(config) < config.run.mut_prob
        #valid_types = filter(t -> t != site.type, [instances(ProteinPropsMod.ProteinType)...])
        site.type = Random.rand(config.rng, instances(ProteinPropsMod.ProteinType))
    end
    if RandUtilsMod.rand_float(config) < config.run.mut_prob
        site.tag = UInt8(RandUtilsMod.rand_int(config, 0, 255))
    end
    if RandUtilsMod.rand_float(config) < config.run.mut_prob
        #valid_actions = filter(t -> t != site.action, [instances(ProteinPropsMod.ProteinAction)...])
        site.action = Random.rand(config.rng, instances(ProteinPropsMod.ProteinAction))
    end
    if RandUtilsMod.rand_float(config) < config.run.mut_prob
        site.arg = RandUtilsMod.rand_int(config, 0, 127)
    end

    #mutate site.threshold and site.consum_rate
    #mutate_floats(config, site, ea_step)
end

function mutate_bind_site(config::Config, site::BindSite, ea_step::Int64)
    #note: all bind sites must not have type Application
    if RandUtilsMod.rand_float(config) < config.run.mut_prob
        valid_types = filter(t -> t != ProteinPropsMod.Application, [instances(ProteinPropsMod.ProteinType)...])
        site.type = Random.rand(config.rng, valid_types)
    end
    if RandUtilsMod.rand_float(config) < config.run.mut_prob
        site.tag = UInt8(RandUtilsMod.rand_int(config, 0, 255))
    end

    #mutate site.threshold and site.consum_rate
    #mutate_floats(config, site, ea_step)
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
    arg::Union{Array{UInt8, 1}, Nothing}=nothing
)
    i = 1
    enum_info = (
        (ProteinPropsMod.ProteinType, type, :type),
        (ProteinPropsMod.ProteinLoc, loc, :loc),
        (ProteinPropsMod.ProteinAction, action, :action)
    )
    while i <= length(enum_info)
        if RandUtilsMod.rand_float(config) < config.run.mut_prob
            enum, options, fieldname = enum_info[i]
            if options == nothing
                #make sure we don't select the current value
                valid_options = filter(val -> val != getfield(props, fieldname), [instances(enum)...])
                new_val = Random.rand(config.rng, valid_options)
            else
                #make sure we don't select the current value
                valid_options = filter(val -> val != getfield(props, fieldname), options)
                new_val = Random.rand(config.rng, valid_options)
            end
            setfield!(props, fieldname, new_val)
        end
        i += 1
    end

    #mutate props.arg
    if RandUtilsMod.rand_float(config) < config.run.mut_prob
        if arg == nothing
            time_factor = 1.0 - ea_step / config.run.ea_steps
            range = Int64(floor(time_factor * config.run.max_protein_arg))
            delta = RandUtilsMod.rand_int(config, 0, range) - range ÷ 2
            #watch out for overflow
            if props.arg + delta > config.run.max_protein_arg
                props.arg = UInt8(config.run.max_protein_arg)
            elseif props.arg + delta < 0
                props.arg = UInt8(0)
            else
                props.arg = UInt8(props.arg + delta)
            end
        else
            props.arg = Random.rand(config.rng, arg)
        end
    end
end

end

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
    #dup_and_swap(indiv, ea_step)
    #dup_and_mut(indiv, ea_step)
    system_level(indiv, ea_step)
end

function dup_and_swap(indiv::Individual, ea_step::Int64)
    gene_index = 1
    first_new_index = -1
    mutated = false
    while gene_index <= length(indiv.genes)
        if (ea_step <= indiv.config.run.gene_dup_gen_limit &&
            indiv.reg_sim_info.produce_count[gene_index] >= indiv.config.run.gene_dup_count_threshold &&
            length(indiv.genes) < indiv.config.run.max_genes)

            gene = indiv.genes[gene_index]
            if RandUtilsMod.rand_float(gene.config) < gene.config.run.dup_mut_prob
                copy = deepcopy(gene)
                #swap the types and tags of the bind sites and prod sites
                for i in 1:length(gene.bind_sites)
                    bsite = gene.bind_sites[i]
                    psite = gene.prod_sites[i]
                    bsite.type, psite.type = psite.type, bsite.type
                    bsite.tag, psite.tag = psite.tag, bsite.tag
                end

                #insert mutated copy after the src gene
                insert_genes(indiv, [copy], gene_index + 1)

                #update the index of first mutation
                if first_new_index == -1 || gene_index + 1 < first_new_index
                    first_new_index = gene_index + 1
                end

                mutated = true #since we've modified the genome, consider it mutated
                gene_index += 1 #skip over the copy
            end
            
        else
            #otherwise, point mutate the gene
            #if indiv.reg_sim_info.produce_count[gene_index] == 0
            mutated = mutated || point_mutate_gene(indiv, gene_index, ea_step)
            #end
        end
        
        gene_index += 1
    end

    #swap locations
    #mutate_location(indiv.config, indiv.genes)

    #update genome_indices
    # foreach(i -> indiv.genes[i].genome_index = i, 1:length(indiv.genes))

    if mutated
        indiv.last_mod = ea_step
        if first_new_index != -1
            update_genome_indices(indiv, first_new_index)
        end
    end
end

function dup_and_mut(indiv::Individual, ea_step::Int64)
    gene_index = 1
    first_new_index = -1
    mutated = false
    while gene_index <= length(indiv.genes)
        if (ea_step <= indiv.config.run.gene_dup_gen_limit &&
            indiv.reg_sim_info.produce_count[gene_index] >= indiv.config.run.gene_dup_count_threshold &&
            length(indiv.genes) <= indiv.config.run.max_genes - 1)

            gene = indiv.genes[gene_index]
            if RandUtilsMod.rand_float(gene.config) < gene.config.run.dup_mut_prob
                copy = deepcopy(gene)
                point_mutate_gene(indiv, gene_index, ea_step)

                #insert mutated copy after the src gene
                insert_genes(indiv, [copy], gene_index + 1)

                #update the index of first mutation
                if first_new_index == -1 || gene_index + 1 < first_new_index
                    first_new_index = gene_index + 1
                end

                mutated = true #since we've modified the genome, consider it mutated
                gene_index += 1 #skip over the copy
            end
            
        else
            #otherwise, point mutate the gene
            #if indiv.reg_sim_info.produce_count[gene_index] == 0
            mutated = mutated || point_mutate_gene(indiv, gene_index, ea_step)
            #end
        end
        
        gene_index += 1
    end

    #swap locations
    #mutate_location(indiv.config, indiv.genes)

    #update genome_indices
    # foreach(i -> indiv.genes[i].genome_index = i, 1:length(indiv.genes))

    if mutated
        indiv.last_mod = ea_step
        if first_new_index != -1
            update_genome_indices(indiv, first_new_index)
        end
    end
end

function system_level(indiv::Individual, ea_step::Int64)
    val, index = findmax(indiv.reg_sim_info.produce_count)
    if (ea_step <= indiv.config.run.gene_sys_level_gen_limit &&
        val >= indiv.config.run.gene_sys_level_threshold &&
        length(indiv.genes) <= indiv.config.run.max_genes - 2)

        if RandUtilsMod.rand_float(indiv.config) < indiv.config.run.sys_level_mut_prob
            #println("System level mutation at ea_step: $(ea_step)")
            temp_to_self_sustaining(indiv, index)
        end
    else
        for i in 1:length(indiv.genes)
            point_mutate_gene(indiv, i, ea_step)
        end
    end
end

function insert_genes(indiv::Individual, new_genes::Array{Gene, 1}, start::Int64)
    for i in 1:length(new_genes)
        gene_index = start + i - 1
        new_gene = new_genes[i]
        
        insert!(indiv.genes, gene_index, new_gene)
    
        insert!(indiv.cell_tree.root.gene_states, gene_index, GeneState(indiv.config.run, new_gene)) #insert a new GeneState into the root cell
        RegSimInfoMod.insert_new_counts(indiv.reg_sim_info, gene_index, 0) #insert new counts for the gene
        IndividualMod.extend_sensors(indiv, gene_index)

        new_num_concs = length(indiv.genes)
        #insert a new conc into the initial proteins, and into the corresponding copy that is in the root cell
        for init_protein in indiv.initial_cell_proteins
            active_protein = ProteinStoreMod.get(indiv.cell_tree.root.proteins, init_protein.props)
            conc = RandUtilsMod.rand_float(indiv.config)
            insert!(init_protein.concs, gene_index, conc)
            
            #this check is needed since there can be duplicate proteins with identical props in cell.initial_cell_proteins
            if length(active_protein.concs) < new_num_concs
                insert!(active_protein.concs, gene_index, conc)
            end
        end
    end
end

function update_genome_indices(indiv::Individual, start_index::Int64)
    for i in start_index:length(indiv.genes)
        indiv.genes[i].genome_index = i
    end
end

function point_mutate_gene(indiv::Individual, gene_index::Int64, ea_step::Int64)
    mutated = false
    gene = indiv.genes[gene_index]
    tree_size = CellTreeMod.size(indiv.cell_tree)
    #bind sites:
    for site in gene.bind_sites
        mutated = mutated || mutate_bind_site(gene.config, site, ea_step, tree_size)
    end

    #prod sites:
    for site in gene.prod_sites
        mutated = mutated || mutate_prod_site(gene.config, site, ea_step)
    end

    mutated
end

function mutate_location(config::Config, genes::Array{Gene, 1})
    if RandUtilsMod.rand_float(config) < config.run.point_mut_prob
        src = RandUtilsMod.rand_int(config, 1, length(genes))
        delta = Random.rand(config.rng, [1, -1])
        dest = src + delta
        if dest < 1
            dest += length(genes)
        elseif dest > length(genes)
            dest = dest % length(genes)
        end

        genes[src].genome_index, genes[dest].genome_index = genes[dest].genome_index, genes[src].genome_index
        genes[src], genes[dest] = genes[dest], genes[src]
    end
end

function mutate_prod_site(config::Config, site::ProdSite, ea_step::Int64)
    mutated = false
    if RandUtilsMod.rand_float(config) < config.run.point_mut_prob
        #valid_types = filter(t -> t != site.type, [instances(ProteinPropsMod.ProteinType)...])
        site.type = Random.rand(config.rng, instances(ProteinPropsMod.ProteinType))
        mutated = true
    end
    if RandUtilsMod.rand_float(config) < config.run.point_mut_prob
        site.tag = UInt8(RandUtilsMod.rand_int(config, 0, 255))
        mutated = true
    end
    if RandUtilsMod.rand_float(config) < config.run.point_mut_prob
        #valid_actions = filter(t -> t != site.action, [instances(ProteinPropsMod.ProteinAction)...])
        site.action = Random.rand(config.rng, instances(ProteinPropsMod.ProteinAction))
        mutated = true
    end
    if RandUtilsMod.rand_float(config) < config.run.point_mut_prob
        site.arg = RandUtilsMod.rand_int(config, 0, 127)
        mutated = true
    end

    #mutate site.threshold and site.consum_rate
    #mutate_floats(config, site, ea_step)

    mutated
end

function mutate_bind_site(config::Config, site::BindSite, ea_step::Int64, tree_size::Int64)
    mutated = false
    #note: all bind sites must not have type Application
    #note: if tree size is 1, there's no point to having neighbour or diffusion bind sites -
    #so the only type left is internal, and no type mutation is possible
    if RandUtilsMod.rand_float(config) < config.run.point_mut_prob && tree_size > 1
        valid_types = filter(t -> t != ProteinPropsMod.Application, [instances(ProteinPropsMod.ProteinType)...])
        site.type = Random.rand(config.rng, valid_types)
        mutated = true
    end
    if RandUtilsMod.rand_float(config) < config.run.point_mut_prob
        site.tag = UInt8(RandUtilsMod.rand_int(config, 0, 255))
        mutated = true
    end

    #mutate site.threshold and site.consum_rate
    #mutate_floats(config, site, ea_step)

    mutated
end

function mutate_floats(config::Config, site::BindSite, ea_step::Int64)
    time_factor = 1.0 - ea_step / config.run.ea_steps
    
    for fieldname in (:threshold, :consum_rate)
        if RandUtilsMod.rand_float(config) < config.run.point_mut_prob
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
        if RandUtilsMod.rand_float(config) < config.run.point_mut_prob
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
    if RandUtilsMod.rand_float(config) < config.run.point_mut_prob
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

#note: right now this will only work with 1 bind site and 1 prod site
function temp_to_self_sustaining(indiv::Individual, index::Int64)
    src_gene = indiv.genes[index]
    src_a_bind_site = src_gene.bind_sites[1]
    src_c_prod_site = src_gene.prod_sites[1]
    #build a list of values that b's tag could take on
    #these include all the values between 1:255 except a's tag and c's tag
    low_exception = min(src_c_prod_site.tag, src_a_bind_site.tag)
    high_exception = max(src_c_prod_site.tag, src_a_bind_site.tag)
    b_tags = vcat(
        collect(1:low_exception - 1),
        collect(low_exception + 1, high_exception - 1),
        collect(high_exception + 1, 2^8-1)
    )
    b_type = ProteinPropsMod.Internal
    b_tag = Random.rand(src_gene.config.rng, b_tags)
    b_action = RandUtilsMod.rand_enum_val(src_gene.config, ProteinPropsMod.ProteinAction)
    b_arg = Random.rand(Int8)

    #note: this does not need to be unique
    a_arg = Random.rand(Int8)
    
    new_gene = Gene(
        src_gene.config,
        src_gene.genome_index + 1, #leave room for the left gene
        GeneMod.Id,
        src_gene.bind_sites,
        Array{ProdSite, 1}([
            ProdSite( #will produce intermediate b protein
                      b_type,
                      b_tag,
                      b_action,
                      b_arg,
                      IndividualMod.initial_threshold,
                      IndividualMod.initial_consum_rate
            )
        ])
    )

    left_gene = Gene(
        src_gene.config,
        src_gene.genome_index,
        GeneMod.Id,
        Array{BindSite, 1}([
            BindSite( #accepts intermediate b protein
                      b_type,
                      b_tag,
                      IndividualMod.initial_threshold,
                      IndividualMod.initial_consum_rate
            )
        ]),
        Array{ProdSite, 1}([
            ProdSite( #produces 'a' protein
                      src_a_bind_site.type,
                      src_a_bind_site.tag,
                      src_a_bind_site.action,
                      a_arg,
                      IndividualMod.initial_threshold,
                      IndividualMod.initial_consum_rate
            )
        ])
    )

    right_gene = Gene(
        src_gene.config,
        src_gene.genome_index + 2,
        GeneMod.Id,
        Array{BindSite, 1}([
            BindSite( #accepts intermediate b protein
                      b_type,
                      b_tag,
                      IndividualMod.initial_threshold,
                      IndividualMod.initial_consum_rate
            )
        ]),
        src_gene.prod_sites
    )

    #left replaces src_gene at index
    indiv.genes[index] = left
    #new gene gets inserted next
    insert!(indiv.genes, index + 1, new_gene)
    #followed by right
    insert!(indiv.genes, index + 2, right)
end

end

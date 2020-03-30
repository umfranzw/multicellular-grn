module MutateMod

using RunMod
using IndividualMod
using GeneMod
using ProteinMod
using ProteinPropsMod
using RandUtilsMod

function mutate(pop::Array{Individual, 1})
    for indiv in pop
        mutate_indiv(indiv)
    end
end

function mutate_indiv(indiv::Individual)
    app_contrib_genes = ChainGraphMod.get_app_contributing_genes(indiv.chain_graph)

    copy_locations = Array{Int64, 1}()
    i = 1
    while i <= length(indiv.genes)
        if indiv.genes[i] in app_contrib_genes
            mut_copy = dup_mutate_gene(indiv.genes[i])
            if mut_copy != nothing
                insert!(indiv.genes, i + 1, mut_copy) #insert mutated copy after the src gene
                push!(copy_locations, i + 1)
                i += 1 #skip over the copy
            end
        else
            point_mutate_gene(indiv.genes[i])
        end
        
        i += 1
    end

    #insert values into the initial protein concs over the spots where the duplicated genes were inserted
    for loc in copy_locations
        for protein in indiv.initial_cell_proteins
            insert!(protein.concs, loc, RandUtilsMod.rand_float(indiv.config))
        end
    end
end

function dup_mutate_gene(gene::Gene)
    copy = nothing
    if RandUtilsMod.rand_float(gene.config) < config.run.mut_prob
        copy = deepcopy(gene)
        point_mutate_gene(copy)
    end

    copy
end

function point_mutate_gene(gene::Gene)
    #Regulatory sites:
    #type & target must stay as they are, other attributes can mutate
    for site in gene.reg_sites
        mutate_props(gene.config, site, [ProteinPropsMod.ProteinRegAction], true)
    end

    #Production sites:
    #target must stay as it is, other attributes can mutate
    for site in gene.prod_sites
        mutate_props(gene.config, site, [ProteinPropsMod.ProteinType, ProteinPropsMod.ProteinRegAction], true)
    end
end

#note: just because a mutation happens doesn't necessarily mean anything has changed - it's possible for random mutations to choose the same val that currently exists...
function mutate_props(
    config::Config,
    props::ProteinProps,
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
            delta = RandUtilsMod.rand_int(config, -config.run.max_arg_mut_step, config.run.max_arg_mut_step)
            props.arg = Int8((props.arg + delta) % 2^7)
        else
            props.arg = Random.rand(config.rng, arg)
        end
        mutated = true
    end

    mutated
end

end

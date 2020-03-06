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
    
    #Commented out because this creates two moving targets that need to be aligned
    # for protein in indiv.initial_cell_proteins
    #     mutate_initial_protein(protein)
    # end
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

function mutate_initial_protein(protein::Protein)
    #type and target must stay the same, others can be mutated
    mutate_props(protein.config, protein.props, [ProteinPropsMod.ProteinRegAction], true)
    
    #concs can also mutate
    for i in 1:length(protein.concs)
        if RandUtilsMod.rand_float(protein.config) < protein.config.run.mut_prob
            r = RandUtilsMod.rand_float(protein.config) #[0, 1]
            r *= (protein.config.run.max_conc_mut) #[0, max_conc_mut]
            r -= protein.config.run.max_conc_mut / 2 #[-max_conc_mut / 2, max_conc_mut / 2]

            protein.concs[i] = clamp(protein.concs[i] + r, 0, 1)
        end
    end
end

#the enums in the given array are the ones that will be (potentially) mutated
#use symbols instead!!!
function mutate_props(config::Config, props::ProteinProps, enums::Array{DataType, 1}, app_action::Bool)
    for enum in enums
        if RandUtilsMod.rand_float(config) < config.run.mut_prob
            new_val = RandUtilsMod.rand_enum_val(config, enum)

            field_types = fieldtypes(ProteinProps)
            field_index = findfirst(e -> e == enum, field_types)
            field_names = fieldnames(ProteinProps)
            name = field_names[field_index]
            setfield!(props, name, new_val)
        end
    end

    if app_action
        props.app_action = UInt8(RandUtilsMod.rand_int(config, 1, Int64(ProteinPropsMod.num_app_actions)))
    end
end

end

module MutateMod

using RunMod
using IndividualMod
using GeneMod
using ProteinMod
using ProteinPropsMod
using RandUtilsMod

function mutate(run::Run, pop::Array{Individual, 1})
    for indiv in pop
        for gene in indiv.genes
            mutate_gene(gene)
        end
        
        for protein in indiv.initial_cell_proteins
            mutate_initial_protein(protein)
        end
    end
end

function mutate_gene(gene::Gene)
    #Regulatory sites:
    #type & target must stay as they are, other attributes can mutate
    for site in gene.reg_sites
        mutate_props(gene.run, site, [ProteinPropsMod.ProteinRegAction, ProteinPropsMod.ProteinAppAction])
    end

    #Production sites:
    #target must stay as it is, other attributes can mutate
    for site in gene.prod_sites
        mutate_props(gene.run, site, [ProteinPropsMod.ProteinType, ProteinPropsMod.ProteinRegAction, ProteinPropsMod.ProteinAppAction])
    end
end

function mutate_initial_protein(protein::Protein)
    #type and target must stay the same, others can be mutated
    mutate_props(protein.run, protein.props, [ProteinPropsMod.ProteinRegAction, ProteinPropsMod.ProteinAppAction])
    
    #concs can also mutate
    for i in 1:length(protein.concs)
        if RandUtilsMod.rand_float(protein.run) < protein.run.mut_prob
            r = RandUtilsMod.rand_float(protein.run) #[0, 1]
            r *= (protein.run.max_conc_mut) #[0, max_conc_mut]
            r -= protein.run.max_conc_mut / 2 #[-max_conc_mut / 2, max_conc_mut / 2]

            protein.concs[i] = clamp(protein.concs[i] + r, 0, 1)
        end
    end
end

#the enums in the given array are the ones that will be (potentially) mutated
function mutate_props(run::Run, props::ProteinProps, enums::Array{DataType, 1})
    for enum in enums
        if RandUtilsMod.rand_float(run) < run.mut_prob
            new_val = RandUtilsMod.rand_enum_val(run, enum)

            field_types = fieldtypes(ProteinProps)
            field_index = findfirst(e -> e == enum, field_types)
            field_names = fieldnames(ProteinProps)
            name = field_names[field_index]
            setfield!(props, name, new_val)
        end
    end
end

end
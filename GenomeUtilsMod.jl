module GenomeUtilsMod

using IndividualMod
using GeneMod
using GeneStateMod
using RegSimInfoMod
using ProteinStoreMod

import RandUtilsMod

function insert_genes(indiv::Individual, new_genes::Array{Gene, 1}, start::Int64)
    half = (1.0 - indiv.config.run.bind_threshold) / 2
    
    for i in 1:length(new_genes)
        gene_index = start + i - 1
        new_gene = new_genes[i]
        
        insert!(indiv.genes, gene_index, new_gene)
    
        insert!(indiv.cell_tree.root.gene_states, gene_index, GeneState(indiv.config.run, new_gene)) #insert a new GeneState into the root cell
        RegSimInfoMod.insert_new_counts(indiv.reg_sim_info, gene_index, 0) #insert new counts for the gene

        new_num_concs = length(indiv.genes)
        #insert a new conc into the initial proteins, and into the corresponding copy that is in the root cell
        for init_protein in indiv.initial_cell_proteins
            active_protein = ProteinStoreMod.get(indiv.cell_tree.root.proteins, init_protein.props)
            conc = RandUtilsMod.rand_floats(indiv.config, indiv.config.run.bind_threshold + half, 1.0)[1]
            insert!(init_protein.concs, gene_index, conc)
            
            #this check is needed since there can be duplicate proteins with identical props in cell.initial_cell_proteins
            if length(active_protein.concs) < new_num_concs
                insert!(active_protein.concs, gene_index, conc)
            end
        end
    end

    update_genome_indices(indiv, start + length(new_genes))
end

function replace_genes(indiv::Individual, new_genes::Array{Gene, 1}, start::Int64)
    for i in 1:length(new_genes)
        gene_index = start + i - 1
        new_gene = new_genes[i]

        #overwrite the existing gene in the indiv's array
        indiv.genes[gene_index] = new_gene
        #the gene state in the root cell contains a reference to the old gene. Update it.
        indiv.cell_tree.root.gene_states[gene_index].gene = new_gene
        #clear any existing bindings
        GeneStateMod.clear_all_bindings(indiv.cell_tree.root.gene_states[gene_index])
    end
end

function update_genome_indices(indiv::Individual, start_index::Int64=1)
    for i in start_index:length(indiv.genes)
        indiv.genes[i].genome_index = i
    end
end

end

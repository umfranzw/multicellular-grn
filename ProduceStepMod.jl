module ProduceStepMod

using IndividualMod
using CellTreeMod
using GeneStateMod
using ProteinStoreMod
using ProteinPropsMod
using ProteinMod
using CellMod

function run_produce(indiv::Individual)
    CellTreeMod.traverse(cell -> run_produce_for_cell(indiv, cell), indiv.cell_tree)
end

function run_produce_for_cell(indiv::Individual, cell::Cell)
    for gene_index in 1:length(cell.gene_states)
        gene_state = cell.gene_states[gene_index]
        gene = indiv.genes[gene_index]
        #note: rates is an array of the form [(prod_site_index, rate), ...]
        rates = GeneStateMod.get_prod_rates(gene_state)

        for (prod_index, rate) in rates
            run_produce_for_site(cell, gene_index, prod_index, rate)
            indiv.reg_sim_info.produce_count[gene_index] += 1
        end
    end
end

function run_produce_for_site(cell::Cell, gene_index::Int64, prod_index::Int64, rate::Float64)
    #get the props for the protein that will be produced
    gene = cell.gene_states[gene_index].gene
    prod_site = gene.prod_sites[prod_index]
    props = ProteinProps(
        prod_site.type,
        prod_site.tag,
        prod_site.action,
        prod_site.arg
    )

    #check if protein already exists in this cell's store
    protein = ProteinStoreMod.get(cell.proteins, props)
    #if not, create and insert it
    if protein == nothing
        if ProteinStoreMod.num_proteins(cell.proteins) < cell.config.run.max_proteins_per_cell
            #note: protein will be initialized with conc values of zero
            #@info @sprintf("Produced protein: %s", props)
            #note: remember to deepcopy the props, so that if this prod site is mutated, the protein's props don't also
            protein = Protein(cell.config, deepcopy(props), false, false, length(cell.gene_states), cell.id)
            ProteinStoreMod.insert(cell.proteins, protein)
        end
    end

    if protein != nothing
        #increment the conc using the rate
        #note: this will only increment the conc directly over the gene
        #the diffusion will spread this out later
        protein.concs[gene_index] = min(protein.concs[gene_index] + rate, 1.0)
    end
end

end

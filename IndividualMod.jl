module IndividualMod

using RunMod
using GeneMod
using GeneStateMod
using ProteinMod
using CellMod
using ProteinStoreMod

import Random
import RandUtilsMod

export Individual,
    rand_init, run_bind

struct Individual
    run::Run
    genes::Array{Gene, 1}
    cells::Array{Cell, 1}
    initial_cell_proteins::Array{Protein, 1}
    protein_store::ProteinStore
end

function rand_init(run::Run)
    genes = map(i -> GeneMod.rand_init(run, i), 1:run.num_genes)
    store = ProteinStore(run)
    initial_cell = Cell(run, genes)
    
    initial_proteins = Array{Protein, 1}()
    seq_vals = Random.shuffle(0:2 ^ ProteinMod.num_bits - 1)
    for val in seq_vals[1:min(length(seq_vals), run.num_initial_proteins)]
        seq = BitArray(val)
        protein = Protein(run, seq, 1, true)
        push!(initial_proteins, protein)

        #also push it in to the store so the first cell has access to it
        #note: we push a copy so the indiv's array stays intact as the simulation modifies proteins
        ProteinStoreMod.insert_protein(store, ProteinMod.copy(protein), nothing)
    end
    
    Individual(run, genes, [initial_cell], initial_proteins, store)
end

function divide_cell(indiv::Individual, src_cell::Cell, dest_index::Int64)
    dest_cell = Cell(indiv.run, indiv.genes)

    ProteinStoreMod.insert_cell(indiv.protein_store, dest_index)
    insert!(indiv.cell, dest_index)

    #insert the initial proteins for the new cell
    #src cell donates half its intra-cell proteins' concentrations to the dest cell
    intra_proteins = ProteinStoreMod.get_proteins_by_scope(indiv.protein_store, ProteinMod.IntraCell)

    for protein in intra_proteins
        protein.concs[src_cell] /= 2
        protein.concs[dest_cell] = copy(protein.concs[src_cell])
    end

    src_cell.energy /= 2
end

function run_growth(indiv::Individual)
    #as we insert, the list will grow - we need to compensate for this
    cells = copy(indiv.cells) #we will loop through this copy so that we don't run into cells inserted on a previous iter
    index_offset = 0 #offset for cells already inserted

    for i in 1:length(cells)
        prob, dir = calc_cell_divide_prob(cells[i])
        if RandUtils.rand_float(indiv.run) < prob
            if dir == ProteinMod.GrowUp
                divide_cell(indiv, cells[i], max(i + index_offset - 1, 1))
            elseif dir == ProteinMod.GrowDown
                divide_cell(indiv, cells[i], min(i + index_offset + 1, length(indiv.cells)))
            end
            index_offset += 1
        end
    end
end

function run_bind(indiv::Individual)
    for i in 1:length(indiv.cells)
        run_bind_for_cell(indiv, i)
    end
end

function run_produce(indiv::Individual)
    for i in 1:length(indiv.cells)
        run_produce_for_cell(indiv, i)
    end
end

function run_regulate(indiv::Individual)
    for i in 1:length(indiv.cells)
        run_regulate_for_cell(indiv, i)
    end
end

function run_diffuse(indiv::Individual)
    DiffusionMod.diffuse_intra_cell_proteins(indiv)
    DiffusionMod.diffuse_inter_cell_proteins(indiv)
end

function run_regulate_for_cell(indiv::Individual, cell_index::Int64)
    cell = indiv.cells[cell_index]
    for gene_index in 1:indiv.run.num_genes
        gene_state = cell.gene_states[gene_index]

        bound_protein = gene_state.reg_site_binding
        if bound_protein != nothing
            affinity = ProteinMod.get_bind_affinity(bound_protein)
            if affinity == ProteinMod.Reg
                reg_action = ProteinMod.get_reg_action(bound_protein)
                if reg_action == ProteinMod.RateUp
                    gene_state.prod_rate = min(gene_state.prod_rate + indiv.run.prod_rate_incr, 1.0)
                elseif reg_action == ProteinMod.RateDown
                    gene_state.prod_rate = max(gene_state.prod_rate - indiv.run.prod_rate_incr, 0.0)
                end
            end
        end
    end
end

function run_produce_for_cell(indiv::Individual, cell_index::Int64)
    cell = indiv.cells[cell_index]
    for gene_index in 1:indiv.run.num_genes
        gene_state = cell.gene_states[gene_index]
        gene = indiv.genes[gene_index]
        
        for site_index in 1:indiv.run.num_bind_sites
            #check if something's bound to the bind site and nothing's bound to the prod site
            bound_protein = gene_state.bind_site_bindings[site_index]
            inhibit_protein = gene_state.prod_site_bindings[site_index]
            if bound_protein != nothing && inhibit_protein == nothing
                product_seq = gene.bind_sites[site_index]
                
                #if product is not yet in store, create it
                product_protein = ProteinStoreMod.get_protein(indiv.protein_store, product_seq)
                if product_protein == nothing
                    product_protein = ProteinMod.Protein(run, product_seq, length(indiv.cells), false)
                    ProteinStoreMod.insert_protein(indiv.protein_store, product_protein, cell)

                #otherwise, make sure this cell is a src
                else
                    #note: don't need to worry about inserting duplicates here, since these are stored in a set
                    ProteinStoreMod.add_src_cell(indiv.protein_store, cell)
                end

                #increment conc above this gene (diffusion will spread it out later)
                product_protein.concs[cell_index][gene_index] = min(produce_protein.concs[cell_index][gene_index] + gene_state.output_rate, 1.0)
            end
        end
    end
end

function run_bind_for_cell(indiv::Individual, cell_index::Int64)
    for gene_index in 1:indiv.run.num_genes
        #regulatory site
        run_binding_for_site(indiv, cell_index, gene_index, indiv.genes[gene_index].reg_site, GeneStateMod.RegSite, indiv.run.reg_bind_threshold)

        #growth site
        run_binding_for_site(indiv, cell_index, gene_index, indiv.genes[gene_index].growth_site, GeneStateMod.GrowthSite, indiv.run.growth_bind_threshold)
        
        for site_index in 1:length(site_seqs)
            #bind sites
            run_binding_for_site(indiv, cell_index, gene_index, indiv.genes[gene_index].bind_sites[site_index], GeneStateMod.BindSite, indiv.run.bind_bind_threshold)
            
            #prod sites
            run_binding_for_site(indiv, cell_index, gene_index, indiv.genes[gene_index].prod_sites[site_index], GeneStateMod.ProdSite, indiv.run.prod_bind_threshold)
        end
    end
    
end

function run_bind_for_site(indiv::Individual, cell_index::Int64, gene_index::Int64, site_seq::BitArray{1}, site_type::GeneStateMod.SiteType, bind_threshold::Float64)
    #build a list of the proteins that can bind to the specified site
    eligigble_proteins = Array{Protein, 1}()
    internal_proteins = values(ProteinStore.get_proteins_by_target(indiv.protein_store, ProteinMod.Internal))
    for protein in internal_proteins
        above_thresh = protein.concs[cell_index][gene_index] >= bind_threshold
        num_diff_bits = ProteinMod.num_bits - BitUtilsMod.count_common_bits(protein.seq, site_seq)
        enough_bit_similarity = num_diff_bits <= run.binding_seq_play

        #we need to make sure that inter-cell proteins don't bind to genes in the cell that produced them (don't "self-bind")
        is_self_binding = ProteinMod.get_scope(protein) == ProteinMod.InterCell && ProteinStore.has_src_cell(indiv.store, protein, indiv.cells[cell_index])

        if above_thresh && enough_bit_similarity && !is_self_binding
            push!(eligable_proteins, protein)
        end
    end

    if length(eligible_proteins) > 0
        conc_sum = foldl((s, p) -> s + p.concs[cell_index][gene_index], eligible_proteins; init=0.0)
        next_bound = 0.0
        wheel = []
        for protein in eligible_proteins
            next_bound = protein.concs[cell_index][gene_index] / conc_sum
            push!(wheel, next_bound)
        end

        wheel[end] = 1.0 #just in case we've got some floating point error

        sel_index = 1
        r = RandUtilsMod.rand_float(run) #random value in [0, 1)
        while r >= wheel[sel_index]
            sel_index += 1
        end

        sel_protein = eligible_proteins[sel_index]

        GeneState.bind(gene_state, sel_protein, site_type, site_index)

    else
        GeneState.unbind(gene_state, site_type, site_index)
    end
end

end

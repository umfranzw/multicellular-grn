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
    rand_init, run_binding

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
        concs = [RandUtilsMod.rand_floats(run, run.num_genes)]
        protein = Protein(run, seq, concs)
        push!(initial_proteins, protein)

        #also push it in to the store so the first cell has access to it
        ProteinStoreMod.insert_protein(store, protein, nothing)
    end
    
    Individual(run, genes, [initial_cell], initial_proteins, store)
end

function run_binding(indiv::Individual)
    for i in 1:length(indiv.cells)
        run_binding(indiv, i)
    end
end

function is_protein_bind_eligible(protein::Protein, cell_index::Int64, gene_index::Int64, site_seq::BitArray{1})
    above_thresh = protein.concs[cell_index][gene_index] >= indiv.run.bind_threshold
    num_diff_bits = ProteinMod.num_bits - BitUtilsMod.count_common_bits(protein.seq, site_seq)
    enough_bit_similarity = num_diff_bits <= run.binding_seq_play

    #we need to make sure that inter-cell proteins don't bind to genes in the cell that produced them (don't "self-bind")
    #self_bind_cond = (ProteinMod.get_scope(protein) == ProteinMod.IntraCell || (cell_index ∉ protein.src_cells))

    above_thresh && enough_bit_similarity# && self_bind_cond
end

function run_binding(indiv::Individual, cell_index::Int64)
    #run binding to regulatory sites
    run_binding(indiv, cell_index, [indiv.genes.reg_site], GeneStateMod.RegSite)
    
    #run binding to bind sites
    run_binding(indiv, cell_index, indiv.genes.bind_sites, GeneStateMod.BindSite)
    
    #run binding to prod sites
    run_binding(indiv, cell_index, indiv.genes.prod_sites, GeneStateMod.ProdSite)
end

function run_binding(indiv::Individual, cell_index::Int64, site_seqs::Array{BitArray{1}, 1}, site_type::GeneStateMod.SiteType)
    for gene_index in 1:indiv.run.num_genes
        for site_index in 1:length(site_seqs)
            eligible_proteins = filter(
                p -> is_protein_bind_eligible(protein, cell_index, gene_index, site_seqs[bind_index]),
                values(ProteinStore.get_proteins_by_target(indiv.protein_store, ProteinMod.Internal))
            )

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
end

end

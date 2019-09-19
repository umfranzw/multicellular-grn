module IndividualMod

using RunMod
using GeneMod
using ProteinMod
using CellMod

import Random
import RandUtilsMod

export Individual,
    rand_init, run_binding

struct Individual
    run::Run
    genes::Array{Gene, 1}
    cells::Array{Cell, 1}
    initial_cell_proteins::Array{Protein, 1}
end

function rand_init(run::Run)
    genes = map(i -> GeneMod.rand_init(run, i), 1:run.num_genes)
    initial_proteins = Array{Protein, 1}()

    seq_vals = Random.shuffle(0:2 ^ ProteinMod.num_bits - 1)
    for val in seq_vals[1:min(length(seq_vals), run.num_initial_proteins)]
        bits = BitArray(val)
        concs = RandUtilsMod.rand_floats(run, run.num_genes)
        protein = Protein(run, bits, -1, concs) #src_gene_index is -1 here to indicate that this is an initial protein (ie. it was not produced by a gene)
        push!(initial_proteins, protein)
    end
    
    initial_cell = Cell(run, genes, initial_proteins)
    
    Individual(run, genes, [initial_cell], initial_proteins)
end

function run_binding(indiv::Individual)
    #intra-cell
    for cell in indiv.cells
        run_intra_cell_binding(indiv, cell)
    end

    #inter-cell
end

function is_protein_bind_eligable(protein::Protein, pos::Int64, site_seq::BitArray{1})
    above_thresh = protein.concs[pos] >= indiv.run.bind_threshold
    num_diff_bits = ProteinMod.num_bits - BitUtilsMod.count_common_bits(protein.seq, site_seq)
    enough_bit_similarity = num_diff_bits <= run.binding_seq_play

    above_thresh && enough_bit_similarity
end

function run_intra_cell_binding(indiv::Individual, cell::Cell)
    for gene_index in 1:indiv.run.num_genes
        for site_index in 1:indiv.run.num_bind_sites
            eligable_proteins = filter(
                p -> is_protein_bind_eligable(p, pos, gene.bind_sites[bind_index]),
                values(cell.proteins[ProteinMod.IntraCell])
            )

            if length(eligable_proteins) > 0
                conc_sum = foldl((s, p) -> s + p.concs[gene_index], eligable_proteins; init=0.0)
                next_bound = 0.0
                wheel = []
                for protein in eligable_proteins
                    next_bound = protein.concs[gene_index] / conc_sum
                    push!(wheel, next_bound)
                end

                wheel[end] = 1.0 #just in case we've got some floating point error

                sel_index = 1
                r = RandUtilsMod.rand_float(run) #random value in [0, 1)
                while r >= wheel[sel_index]
                    sel_index += 1
                end

                sel_protein = eligable_proteins[sel_index]

                GeneState.bind(gene_state, sel_protein, site_index)

            else
                GeneState.unbind(gene_state, site_index)
            end
            
        end
    end
    
end

end

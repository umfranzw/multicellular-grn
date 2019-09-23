module DiffusionMod

using IndividualMod
using ProteinStoreMod
using ProteinMod

function get(array, row, col, num_genes, num_cells)
    if row < 1 || col < 1 || row > height || col > width
        return 0.0
    else
        return array[row][col]
    end
end

function diffuse_intra_cell_proteins(indiv::Individual)
    cols = indiv.run.num_genes
    intra_cell_proteins = ProteinStoreMod.get_proteins_by_scope(indiv.protein_store, ProteinMod.IntraCell)
    
    for protein in intra_cell_proteins
        row_set = ProteinStoreMod.get_src_cells(indiv.protein_store, protein)
        for i in row_set
            new_row = zeros(Float64, cols)
            for j in 1:cols
                new_row[j] = (1 - 4 * indiv.run.diff_dt * indiv.run.diff_alpha / indiv.run.diff_h^2) * protein.concs[i][j] + indiv.run.diff_dt * indiv.run.diff_alpha * ((get(protein.concs, i, j - 1) + get(protein.concs, i - 1, j) + get(protein.concs, i + 1, j) + get(protein.concs, i, j + 1)) / indiv.run.diff_h^2)
            end
            protein.concs[i] = new_row
        end
    end
end

function diffuse_inter_cell_proteins(indiv::Individual)
    rows = length(indiv.cells)
    cols = indiv.run.num_genes
    inter_cell_proteins = ProteinStoreMod.get_proteins_by_scope(indiv.protein_store, ProteinMod.InterCell)

    for protein in inter_cell_proteins
        new_concs = map(i -> zeros(Float64, rows), 1:num_rows)
        for i in 1:rows
            for j in 1:cols
                new_concs[i][j] = (1 - 4 * indiv.run.dt * indiv.run.alpha / indiv.run.h^2) * protein.concs[i][j] + indiv.run.dt * indiv.run.alpha * ((get(protein.concs, i, j - 1) + get(protein.concs, i - 1, j) + get(protein.concs, i + 1, j) + get(protein.concs, u, j + 1)) / indiv.run.h^2)
            end
        end
        protein.concs = new_concs
    end
end

end

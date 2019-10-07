module DiffusionMod

using IndividualMod
using ProteinStoreMod
using ProteinMod
using CellTreeMod

#note: row_offset must be in [-1, 0, 1]
function get_2D(cell::Cell, protein::Protein, row_offset::Int64, col::Int64, width::Int64)
    outside = col < 1 || col > width
    outside = outside || (row_offset == -1 && cell.parent == nothing)
    outside = outside || (row_offset == 1 && length(cell.children) == 0)

    if outside
        return 0.0
    else
        if row_offset == -1
            num_siblings = length(cell.parent.children)
            return get_from_cell(cell.parent, protein.props, col) / num_siblings
            
        elseif row_offset == 0
            return protein.concs[col]
            
        else
            sum = 0.0
            for child in cell.children
                sum += get_from_cell(child, protein.props, col)
            end

            return sum
        end
    end
end

function get_from_cell(cell::Cell, props::ProteinProps, col::Int64)
    protein = ProteinStoreMod.get(cell.protein_store, props)
    if protein == nothing
        return 0.0
    else
        return protein.concs[col]
    end
end

function get_1D(array::Array{Float64, 1}, col::Int64, len::Int64)
    if col < 1 || col > len
        return 0.0
    else
        return array[col]
    end
end

function diffuse_intra_cell_proteins(indiv::Individual)
    CellTreeMod.traverse(indiv.root_cell, cell -> diffuse_intra_cell_proteins(indiv, cell))
end

function diffuse_intra_cell_proteins(indiv::Individual, cell::Cell)
    cols = indiv.run.num_genes
    intra_cell_proteins = ProteinStoreMod.get_proteins_by_target(cell.protein_store, ProteinMod.Intra)
    
    for protein in intra_cell_proteins
        new_concs = zeros(Float64, cols)
        
        for j in 1:cols
            new_concs[j] = (1 - 4 * indiv.run.diff_dt * indiv.run.diff_alpha / indiv.run.diff_h^2) * protein.concs[j] + indiv.run.diff_dt * indiv.run.diff_alpha * ((get_1D(protein.concs, j - 1, cols) + get_1D(protein.concs, j, cols) + get_1D(protein.concs, j, cols) + get_1D(protein.concs, j + 1, cols)) / indiv.run.diff_h^2)
        end
        protein.concs = new_concs
    end
end

function diffuse_inter_cell_proteins(indiv::Individual)
    cols = indiv.run.num_genes
    inter_cell_proteins = ProteinStoreMod.get_proteins_by_target(cell.protein_store, ProteinMod.Inter)
    
    for protein in inter_cell_proteins
        results_dict = Dict{Cell, Array{Float64, 1}}()
        #this will recursively traverse the cell tree, putting the concs for the next timestep into results_dict
        diffuse_inter_cell_proteins(indiv.root_cell, protein.props, cols, results_dict)

        #update the current protein's concs in each cell in the tree
        CellTreeMod.traverse(cell, c -> update_concs(c, results_dict))
    end
end

function diffuse_inter_cell_proteins(cell::Cell, props::ProteinProps, cols::Int64, results_dict::Dict{Cell, Array{Float64, 1}})
    protein = ProteinStoreMod.get(cell.protein_store, props)
    #if the protein that's diffusing doesn't exist in this cell, create it (with zeroed concs).
    #Note: If the diffusion results in very low concs, the decay step will remove it later.
    if protein == nothing
        protein = Protein(cell.run, props, false)
        ProteinStoreMod.insert(cell.protein_store, protein, false)
    end
    
    new_concs = zeros(Float64, cols)
    for j in 1:cols
        new_concs[j] = (1 - 4 * indiv.run.dt * indiv.run.alpha / indiv.run.h^2) * protein.concs[j] +
            indiv.run.dt * indiv.run.alpha * ((get_2D(cell, protein, 0, j - 1, width) +
                                               get_2D(cell, protein, -1, j, width) +
                                               get_2D(cell, protein, 1, j, width) +
                                               get_2D(cell, protein, 0, j + 1, width)) /
                                              indiv.run.h^2)
        #make sure the value stays in [0.0, 1.0]
        #note: not sure if this is guarenteed here because of the way we're diffusing between parents and children in the tree structure (it would be in a matrix)...
        new_concs[j] = clamp(new_concs[j], 0.0, 1.0)
    end

    #note: we cannot update protein.concs until the whole tree has been traversed (otherwise subsequent calls may use the updated value instead of the original), since this cell may be referred to more than once. So we save it in this dictionary instead.
    results_dict[cell] = new_concs

    for child in cell.children
        diffuse_inter_cell_proteins(child, props, cols, results_Dict)
    end
end

function update_concs(cell::Cell, props::ProteinProps, results_dict::Dict{Cell, Array{Float64, 1}})
    #note: the protein should always be present here, since we inserted it if it wasn't there in diffuse_inter_cell_proteins()
    protein = ProteinStoreMod.get(cell.protein_store, props)
    protein.concs = results_dict[cell]
end

end

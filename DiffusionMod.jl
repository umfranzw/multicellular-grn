module DiffusionMod

using RunMod
using ProteinStoreMod
using ProteinMod
using ProteinPropsMod
using CellTreeMod
using CellMod
using Printf

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
    protein = ProteinStoreMod.get(cell.proteins, props)
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

function diffuse_intra_cell_proteins(cell_tree::CellTree)
    CellTreeMod.traverse(cell -> diffuse_intra_cell_proteins_for_cell(cell), cell_tree)
end

function diffuse_intra_cell_proteins_for_cell(cell::Cell)
    cols = cell.run.num_genes
    intra_cell_proteins = ProteinStoreMod.get_by_target(cell.proteins, ProteinPropsMod.Intra)
    
    for protein in intra_cell_proteins
        new_concs = zeros(Float64, cols)
        
        for j in 1:cols
            new_concs[j] = (1 - 4 * cell.run.diff_dt * cell.run.diff_alpha / cell.run.diff_h^2) * protein.concs[j] + cell.run.diff_dt * cell.run.diff_alpha * ((get_1D(protein.concs, j - 1, cols) + get_1D(protein.concs, j, cols) + get_1D(protein.concs, j, cols) + get_1D(protein.concs, j + 1, cols)) / cell.run.diff_h^2)
        end
        protein.concs = new_concs
    end
end

function get_all_inter_cell_props(cell_tree::CellTree)
    #go through the whole cell tree and build a set of all the inter-cell protein's props
    inter_cell_props = Set{ProteinProps}()
    CellTreeMod.traverse(
        cell -> map(
            p -> push!(inter_cell_props, p.props),
            ProteinStoreMod.get_by_target(cell.proteins, ProteinPropsMod.Inter)
        ),
        cell_tree
    )
    
    inter_cell_props
end

function diffuse_inter_cell_proteins(cell_tree::CellTree)
    #build a set of all of the inter-cell proteins' props (across all cells in the tree)
    inter_cell_props = get_all_inter_cell_props(cell_tree)

    #for each protein, traverse the tree and diffuse it in each cell - storing the new concs in a dictionary
    results = Dict{Cell, Dict{ProteinProps, Array{Float64, 1}}}()

    for props in inter_cell_props
        CellTreeMod.traverse_bf(cell -> diffuse_inter_cell_proteins_for_props(cell, props, results), cell_tree)
    end

    #update the current protein's concs in each cell in the tree
    if length(results) > 0 #note: it's possible there may not be any inter-cellular proteins...
        CellTreeMod.traverse(cell -> update_concs(cell, results), cell_tree)
    end
end

function diffuse_inter_cell_proteins_for_props(cell::Cell, props::ProteinProps, results::Dict{Cell, Dict{ProteinProps, Array{Float64, 1}}})
    cols = cell.run.num_genes
    protein = ProteinStoreMod.get(cell.proteins, props)
    #if the protein that's diffusing doesn't exist in this cell, create it (with zeroed concs).
    #Note: If the diffusion results in zeroed or very low concs, the decay step will remove it later.
    if protein == nothing
        protein = Protein(cell.run, ProteinPropsMod.copy(props), false)
        ProteinStoreMod.insert(cell.proteins, protein, false)
    end
    
    new_concs = zeros(Float64, cols)
    for j in 1:cols
        new_concs[j] = (1 - 4 * cell.run.diff_dt * cell.run.diff_alpha / cell.run.diff_h^2) * protein.concs[j] +
            cell.run.diff_dt * cell.run.diff_alpha * ((get_2D(cell, protein, 0, j - 1, cols) +
                                                       get_2D(cell, protein, -1, j, cols) +
                                                       get_2D(cell, protein, 1, j, cols) +
                                                       get_2D(cell, protein, 0, j + 1, cols)) /
                                                      cell.run.diff_h^2)
        #make sure the value stays in [0.0, 1.0]
        #note: not sure if this is guarenteed here because of the way we're diffusing between parents and children in the tree structure (it would be in a matrix)...
        new_concs[j] = clamp(new_concs[j], 0.0, 1.0)
    end

    #note: we cannot update protein.concs until the whole tree has been traversed (otherwise subsequent
    #calls may use the updated value instead of the original), since this cell may be referred to more than once.
    #So we save it in this dictionary instead.
    if cell ∉ keys(results)
        results[cell] = Dict{ProteinProps, Array{Float64, 1}}()
    end
    results[cell][props] = new_concs
end

function update_concs(cell::Cell, results::Dict{Cell, Dict{ProteinProps, Array{Float64, 1}}})
    #note: the protein should always be present here, since we inserted it if it wasn't there in diffuse_inter_cell_proteins()
    for (props, new_concs) in results[cell]
        protein = ProteinStoreMod.get(cell.proteins, props)
        protein.concs = new_concs
    end
end

end

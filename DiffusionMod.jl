module DiffusionMod

using RunMod
using ProteinStoreMod
using ProteinMod
using ProteinPropsMod
using CellTreeMod
using CellMod
using Printf

struct TreeInfo
    cell_to_level::Dict{Cell, Int64}
    level_to_cell::Dict{Int64, Array{Cell, 1}}

    function TreeInfo(tree::CellTree)
        cell_to_level = Dict{Cell, Int64}()
        level_to_cell = Dict{Int64, Array{Cell, 1}}()
        if tree.root != nothing
            build_info(tree.root, cell_to_level, level_to_cell, 0)
        end

        TreeInfo(cell_to_level, level_to_cell)
    end
end

function build_info(cell::Cell, cell_to_level::Dict{Cell, Int64}, level_to_cell::Dict{Int64, Array{Cell, 1}}, level::Int64)
    cell_to_level[cell] = level
    if level ∉ keys(level_to_cell)
        level_to_cell[level] = Array{Cell}()
    end
    push!(level_to_cell, cell)

    for child in cell.children
        build_info(child, cell_to_level, level_to_cell, level + 1)
    end
end

function diffuse_intra_cell_proteins(cell_tree::CellTree)
    CellTreeMod.traverse(cell -> diffuse_intra_cell_proteins_for_cell(cell), cell_tree)
end

function diffuse_intra_cell_proteins_for_cell(cell::Cell)
    cols = length(cell.gene_states)
    intra_cell_proteins = ProteinStoreMod.get_by_target(cell.proteins, ProteinPropsMod.Intra)
    
    for protein in intra_cell_proteins
        new_concs = zeros(Float64, cols)
        
        for j in 1:cols
            new_concs[j] = (1 - 2 * cell.config.run.diff_dt * cell.config.run.diff_alpha / cell.config.run.diff_h^2) * protein.concs[j] +
                cell.config.run.diff_dt * cell.config.run.diff_alpha * ((get_1D(protein.concs, j - 1, cols) +
                                                                         get_1D(protein.concs, j + 1, cols)) / cell.config.run.diff_h^2)
        end
        protein.concs = new_concs
    end
end

function get_1D(array::Array{Float64, 1}, col::Int64, len::Int64)
    if col < 1 || col > len
        return 0.0
    else
        return array[col]
    end
end

function diffuse_inter_cell_proteins(tree::CellTree)
    #get a set of all of the inter-cell proteins' props (across all cells in the tree)
    inter_proteins = values(tree.inter_proteins)
    info = TreeInfo(tree)

    #for each protein, traverse the tree and diffuse it in each cell - storing the new concs in a dictionary
    results = Dict{ProteinProps, Dict{Cell, Float64}}()

    for protein in inter_proteins
        CellTreeMod.traverse_bf(cell -> diffuse_inter_cell_protein(cell, protein, results, info), tree)
    end

    #update the current protein's concs in each cell in the tree
    if length(results) > 0 #note: it's possible there may not be any inter-cellular proteins...
        CellTreeMod.traverse(cell -> update_concs(cell, results), tree)
    end
end

function diffuse_inter_cell_protein(cell::Cell, protein::InterProtein, results::Dict{ProteinProps, Dict{Cell, Float64}}, info::TreeInfo)
    level = info.cell_to_level[cell]
    num_cols = length(info.level_to_cell[level])
    
    #if the protein that's diffusing doesn't exist in this cell, create it (with zeroed concs).
    #Note: If the diffusion results in zeroed or very low concs, the decay step will remove it later.
    # if protein == nothing
    #     protein = Protein(cell.config, deepcopy(props), false, false, length(cell.gene_states))
    #     ProteinStoreMod.insert(cell.proteins, protein, false)
    # end
    
    new_concs = zeros(Float64, num_cols)
    for j in 1:num_cols
        new_concs[j] = (1 - 4 * cell.config.run.diff_dt * cell.config.run.diff_alpha / cell.config.run.diff_h^2) * protein.concs[j] +
            cell.config.run.diff_dt * cell.config.run.diff_alpha * ((get_2D(cell, protein, 0, j - 1, num_cols) +
                                                       get_2D(cell, protein, -1, j, num_cols) +
                                                       get_2D(cell, protein, 1, j, num_cols) +
                                                       get_2D(cell, protein, 0, j + 1, num_cols)) /
                                                      cell.config.run.diff_h^2)
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

#note: row_offset must be in [-1, 0, 1]
function get_2D(cell::Cell, protein::InterProtein, row_offset::Int64, col::Int64, width::Int64, info::TreeInfo)
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

function update_concs(cell::Cell, results::Dict{Cell, Dict{ProteinProps, Array{Float64, 1}}})
    #note: the protein should always be present here, since we inserted it if it wasn't there in diffuse_inter_cell_proteins()
    for (props, new_concs) in results[cell]
        protein = ProteinStoreMod.get(cell.proteins, props)
        protein.concs = new_concs
    end
end

end

module DiffusionMod

using RunMod
using ProteinStoreMod
using ProteinMod
using ProteinPropsMod
using CellTreeMod
using CellMod
using Printf

function diffuse_intra_cell_proteins(tree::CellTree)
    CellTreeMod.traverse(cell -> diffuse_intra_cell_proteins_for_cell(cell), tree)
end

function get_all_but_diffusion(ps::ProteinStore)
    proteins = Array{Protein, 1}()
    for type in instances(ProteinPropsMod.ProteinType)
        if type != ProteinPropsMod.Diffusion
            append!(proteins, ProteinStoreMod.get_by_type(ps, type))
        end
    end

    proteins
end

function diffuse_intra_cell_proteins_for_cell(cell::Cell)
    cols = length(cell.gene_states)
    intra_cell_proteins = get_all_but_diffusion(cell.proteins)
    
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

function get_all_diffusion_proteins(tree::CellTree)
    proteins = Set{Protein}()
    CellTreeMod.traverse(cell -> append_diffusion_proteins_from_cell(proteins, cell), tree)

    proteins
end

function append_diffusion_proteins_from_cell(set::Set{Protein}, cell::Cell)
    diff_proteins = ProteinStoreMod.get_by_type(cell.proteins, ProteinPropsMod.Diffusion)
    for protein in diff_proteins
        push!(set, protein)
    end
end

function diffuse_inter_cell_proteins(tree::CellTree)
    #get a set of all of the inter-cell proteins' props (across all cells in the tree)
    inter_proteins = get_all_diffusion_proteins(tree)
    info = TreeInfo(tree)

    #for each protein, traverse the tree and diffuse it in each cell - storing the new concs in a dictionary
    for protein in inter_proteins
        results = Dict{Tuple{Int64, Int64}, Array{Float64, 1}}()
        CellTreeMod.traverse_bf(cell -> diffuse_inter_cell_protein(cell, protein, results, info), tree)
        #update the concs
        for ((level, col), new_concs) in results
            #println(new_concs)
            cell = info.level_to_cell[level][col]
            cell_protein = ProteinStoreMod.get(cell.proteins, protein.props)
            if cell_protein != nothing #may not be present if cell already contained max # of proteins and we couldn't insert it
                cell_protein.concs = new_concs
            end
        end
    end
end

function diffuse_inter_cell_protein(cell::Cell, diffusing_protein::Protein, results::Dict{Tuple{Int64, Int64}, Array{Float64, 1}}, info::TreeInfo)
    #if the protein that's diffusing doesn't exist in this cell, create it (with zeroed concs).
    #Note: If the diffusion results in zeroed or very low concs, the decay step will remove it later.
    protein = ProteinStoreMod.get(cell.proteins, diffusing_protein.props)
    if protein == nothing
        if ProteinStoreMod.num_proteins(cell.proteins) < cell.config.run.max_proteins_per_cell
            protein = Protein(cell.config, deepcopy(diffusing_protein.props), false, false, length(cell.gene_states), cell.id)
            ProteinStoreMod.insert(cell.proteins, protein)
        end
    end
    
    if protein != nothing
        cell_level = info.cell_to_level[cell]
        cell_col = findall(n -> n == cell, info.level_to_cell[cell_level])[1]

        new_concs = Array{Float64, 1}()
        num_concs = length(cell.gene_states)
        
        for i in 1:num_concs
            new_conc = (1 - 4 * cell.config.run.diff_dt * cell.config.run.diff_alpha / cell.config.run.diff_h^2) *
                protein.concs[i] +
                cell.config.run.diff_dt * cell.config.run.diff_alpha * ((
                    get_2D(protein, cell_level, cell_col, i, 0, -1, info, num_concs) +
                    get_2D(protein, cell_level, cell_col, i, -1, 0, info, num_concs) +
                    get_2D(protein, cell_level, cell_col, i, 1, 0, info, num_concs) +
                    get_2D(protein, cell_level, cell_col, i, 0, 1, info, num_concs)
                ) / cell.config.run.diff_h^2)
            
            #make sure the value stays in [0.0, 1.0]
            #note: not sure if this is guarenteed here because of the way we're diffusing between parents and children in the tree structure (it would be in a matrix)...
            new_conc = clamp(new_conc, 0.0, 1.0)
            push!(new_concs, new_conc)
        end

        results[(cell_level, cell_col)] = new_concs
    end
end

function get_protein_conc_or_zero(cell::Cell, props::ProteinProps, conc_col::Int64)
    protein = ProteinStoreMod.get(cell.proteins, props)
    if protein == nothing
        conc = 0
    else
        conc = protein.concs[conc_col]
    end

    conc
end

#note: this function will not handle diagonals (i.e. one of {level_offset, col_offset} must be set to 0)
function get_2D(protein::Protein, cell_level::Int64, cell_col::Int64, conc_col::Int64, row_offset::Int64, col_offset::Int64, info::TreeInfo, num_conc_cols::Int64)
    cur_cell = info.level_to_cell[cell_level][cell_col]
    result = 0.0

    #get from above
    if row_offset == -1 && col_offset == 0
        if cur_cell.parent != nothing
            target_conc_col = conc_col + col_offset
            if target_conc_col != 0 && target_conc_col != num_conc_cols + 1 #if not falling off left or right edge
                num_siblings = length(cur_cell.parent.children)
                result = get_protein_conc_or_zero(cur_cell.parent, protein.props, target_conc_col) / num_siblings
            end
        end

    #get from left
    elseif row_offset == 0 && col_offset == -1
        target_conc_col = conc_col + col_offset
        if target_conc_col == 0 #fall off the left edge
            #check to see if we have a left sibling we can grab from
            if cell_col > 1 #have left sibling
                sibling = info.level_to_cell[cell_level][cell_col - 1]
                result = get_protein_conc_or_zero(sibling, protein.props, num_conc_cols) #return the rightmost conc
            end
        else #within bounds in current cell
            result = get_protein_conc_or_zero(cur_cell, protein.props, target_conc_col)
        end

    #get from right
    elseif row_offset == 0 && col_offset == 1
        target_conc_col = conc_col + col_offset
        if target_conc_col == num_conc_cols + 1 #fall off the right edge
            #check to see if we have a right sibling we can grab from
            if cell_col < length(info.level_to_cell[cell_level]) #have right sibling
                sibling = info.level_to_cell[cell_level][cell_col + 1]
                result = get_protein_conc_or_zero(sibling, protein.props, 1) #return the leftmost conc
            end
        else #within bounds
            result = get_protein_conc_or_zero(cur_cell, protein.props, target_conc_col)
        end

    #get from below
    elseif row_offset == 1 && col_offset == 0
        if length(cur_cell.children) > 0
            target_conc_col = conc_col + col_offset
            if target_conc_col != 0 && target_conc_col != num_conc_cols + 1 #if not falling off left or right edge
                sum = 0.0
                for child in cur_cell.children
                    sum += get_protein_conc_or_zero(child, protein.props, target_conc_col)
                end
                result = sum
            end
        end

    else
        error("Diffusion error in get_2D(): Invalid offsets")
    end

    result
end

end

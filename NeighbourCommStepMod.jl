module NeighbourCommStepMod

using IndividualMod
using CellTreeMod
using CellMod
using ProteinStoreMod
using ProteinPropsMod
using ProteinMod

function run_neighbour_comm(indiv::Individual)
    info = TreeInfo(indiv.cell_tree)
    CellTreeMod.traverse(cell -> run_neighbour_comm_for_cell(cell, info), indiv.cell_tree)
end

function get_neighbour_at(cell::Cell, info::TreeInfo, loc::Int64)
    neighbour = nothing

    if loc == Int64(ProteinPropsMod.Top)
        if cell.parent != nothing
            neighbour = cell.parent
        end

    elseif loc == Int64(ProteinPropsMod.Left) || loc == Int64(ProteinPropsMod.Right)
        level = info.cell_to_level[cell]
        row = info.level_to_cell[level]
        cell_index = findall(c -> c == cell, row)[1]

        if loc == Int64(ProteinPropsMod.Left)
            if cell_index > 1
                neighbour = row[cell_index - 1]
            end
        else
            if cell_index < length(row)
                neighbour = row[cell_index + 1]
            end
        end

        #bottom (child)
    else
        num_children = length(cell.children)
        offset = loc - 3 + 1 #- 3 for left, top, right
        if num_children > 0 && offset <= num_children
            neighbour = cell.children[offset]
        end
    end

    neighbour
end

#returns cell's loc with respect to neighbour
function get_loc_relationship(cell::Cell, neighbour::Cell, neighbour_loc::Int64)
    rel = nothing

    if neighbour_loc == Int64(ProteinPropsMod.Top)
        offset = 1
        while offset <= length(neighbour.children) && cell != neighbour.children[offset]
            offset += 1
        end
        if offset <= length(neighbour.children)
            rel = Int8(2 + offset) #left=0, top=1, right=2, first_child=3, ...
        end

    elseif neighbour_loc == Int64(ProteinPropsMod.Left)
        rel = Int8(ProteinPropsMod.Right)

    elseif neighbour_loc == Int64(ProteinPropsMod.Right)
        rel = Int8(ProteinPropsMod.Left)

    else #bottom
        rel = Int8(ProteinPropsMod.Top)
    end

    rel
end

#transfers proteins to cell from neighbours
function run_neighbour_comm_for_cell(cell::Cell, info::TreeInfo)
    for src_loc in 0 : 2 + cell.config.run.max_children #the loc we're taking the proteins from
        neighbour = get_neighbour_at(cell, info, src_loc) #the cell we're taking the proteins from

        if neighbour != nothing
            #get neighbour's neighbour proteins with opposite locs (the ones that are being sent to this cell)
            #note: we use the id to ensure that we only get the ones that originated in the neighbour cell (it's possible they may have been transferred in from another neighbour)
            protein_loc = get_loc_relationship(cell, neighbour, src_loc)
            neighbour_proteins = ProteinStoreMod.get_neighbour_proteins_by_loc(neighbour.proteins, protein_loc, neighbour.id)

            if length(neighbour_proteins) > 0
                #compute the max amount of each protein that we can accept
                #accept_amount = cell.sensors[src_loc] / length(neighbour_proteins)
                for neighbour_protein in neighbour_proteins
                    #transfer_amount = min.(neighbour_protein.concs, accept_amount)
                    transfer_amount = neighbour_protein.concs

                    #add the neighbour protein into the source cell
                    dest_protein = ProteinStoreMod.get(cell.proteins, neighbour_protein.props)
                    if dest_protein == nothing
                        if ProteinStoreMod.num_proteins(cell.proteins) < cell.config.run.max_proteins_per_cell
                            dest_protein = Protein(cell.config, deepcopy(neighbour_protein.props), false, false, length(cell.gene_states), neighbour.id)
                            ProteinStoreMod.insert(cell.proteins, dest_protein)
                        end
                    end
                    if dest_protein != nothing
                        #cell.sensors[src_loc] -= transfer_amount
                        dest_protein.concs = clamp.(dest_protein.concs + transfer_amount, 0.0, 1.0)
                        neighbour_protein.concs -= transfer_amount
                    end

                    #note: if the neighbour_protein was completely tranferred, the decay step will delete it later
                end
            end
        end
    end
end

end

module SearchMod

using DataMod
using CellTreeMod

#Iterates through all cells in the given ranges, calling cell_fcn for each one.
#cell_fcn should accept a cell and return a tuple of the form (Bool, Any), where the Bool
#indicates whether the cell should be included in the results, and the Any is custom user data.
#This function collects the results from the cells for which the function returned True, and returns
#and array of the form [((ea_step, pop_index, reg_step),
#                   (cell_row, cell_col), #in tree
#                   custom_data), #the Any from the cell_fcn return value
#                  ...]
function search_cells(data::Data, ea_step_range::StepRange{Int64, Int64}, pop_index_range::UnitRange{Int64}, reg_step_range::UnitRange{Int64}, cell_fcn::Function)
    results = Array{Tuple{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64}, Any}, 1}()
    
    for ea_step in ea_step_range
        for pop_index in pop_index_range
            for reg_step in reg_step_range
                tree = DataMod.get_tree(data, ea_step, pop_index, reg_step)
                info = TreeInfo(tree)
                levels = sort(collect(keys(info.level_to_cell)))
                for level in levels
                    row_cells = info.level_to_cell[level]
                    for col in 1:length(row_cells)
                        accept, custom_data = cell_fcn(row_cells[col])
                        if accept
                            result = ((ea_step, pop_index, reg_step),
                                      (level, col),
                                      custom_data
                                      )
                            push!(results, result)
                        end
                    end
                end
            end
        end
    end

    results
end

function search_indivs(data::Data, ea_step_range::StepRange{Int64, Int64}, pop_index_range::UnitRange{Int64}, reg_step_range::UnitRange{Int64}, indiv_fcn::Function)
    results = Array{Tuple{Tuple{Int64, Int64, Int64}, Any}, 1}()
    
    for ea_step in ea_step_range
        for pop_index in pop_index_range
            for reg_step in reg_step_range
                indiv = DataMod.get_indiv(data, ea_step, pop_index, reg_step)
                accept, custom_data = indiv_fcn(indiv)
                if accept
                    result = ((ea_step, pop_index, reg_step),
                              custom_data
                              )
                    push!(results, result)
                end
            end
        end
    end

    results
end

function reduce_indivs(data::Data, ea_step_range::StepRange{Int64, Int64}, pop_index_range::UnitRange{Int64}, reg_step_range::UnitRange{Int64}, indiv_fcn::Function)
    last_result = nothing
    
    for ea_step in ea_step_range
        for pop_index in pop_index_range
            for reg_step in reg_step_range
                indiv = DataMod.get_indiv(data, ea_step, pop_index, reg_step)
                last_result = indiv_fcn(indiv, (ea_step, pop_index, reg_step), last_result)
            end
        end
    end

    last_result
end

end

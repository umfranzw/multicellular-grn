using SearchMod
using DataMod
using CellMod
using ProteinPropsMod

function save_cell_results(filename::String, results::Array{Tuple{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64}, Any}, 1}, user_data_fcn::Function)
    file = open(filename, "w")
    write(file, "ea_step, pop_index, reg_step, cell_level, cell_col, user_data\n")
    for row in results
        ((ea_step, pop_index, reg_step), (cell_level, cell_col), user_data) = row
        user_data_str = user_data_fcn(user_data)
        write(file, "$(ea_step), $(pop_index), $(reg_step), $(cell_level), $(cell_col), $(user_data_str)\n")
    end
    
    close(file)
end

function print_cell_results(results::Array{Tuple{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64}, Any}, 1})
    for row in results
        println(row)
    end
end

function cell_with_binding(cell::Cell)
    bind_strs = Array{String, 1}()
    for gs in cell.gene_states
        for binding in gs.bindings
            if binding != nothing
                push!(bind_strs, ProteinPropsMod.to_str(binding.props))
            end
        end
    end

    if length(bind_strs) > 0
        result = (true, bind_strs)
    else
        result = (false, nothing)
    end

    result
end

data = Data("data")
run = DataMod.get_run(data)

ea_step_range = run.step_range
pop_index_range = 1:run.pop_size
reg_step_range = 1:run.reg_steps

results = SearchMod.search_cells(data, ea_step_range, pop_index_range, reg_step_range, cell_with_binding)

DataMod.close(data)
save_cell_results("results.csv", results, bind_strs -> join(bind_strs, ","))

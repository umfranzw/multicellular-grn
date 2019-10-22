module FitnessMod

using IndividualMod
using CellTreeMod

function eval(indiv::Individual)
    expr_str = "f(x) = "
    expr_str *= CellTreeMod.to_expr_str(indiv.cell_tree)

    #f(x) = x * 2
    test_data = [(1, 2), (2, 4)]
    num_passed = 0
    
    for (input, output) in test_data
        try
            test_str = "$(expr_str); f($(input))"
            test_expr = Meta.parse(test_str)
            result = Meta.eval(test_expr)
            num_passed += Int64(result == output)
        catch
            continue
        end
    end
    
    indiv.fitness = 1.0 - num_passed / length(test_data)
end

end

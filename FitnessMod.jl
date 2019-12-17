module FitnessMod

using IndividualMod
using CellTreeMod
using Printf

function eval(indiv::Individual)
    expr_str = "f(x) = "
    expr_str *= CellTreeMod.to_expr_str(indiv.cell_tree)

    #f(x) = x * 2
    #test_data = [(1, 2), (2, 4)]

    #f(x) = x * 2 + 1
    test_data = [(1, 3), (2, 5), (3, 7)]
    
    fitness = 1.0
    chunk = 1 / length(test_data)
    
    for (input, output) in test_data
        try
            test_str = "$(expr_str); f($(input))"
            test_expr = Meta.parse(test_str)
            result = Meta.eval(test_expr)

            if result == output
                fitness -= chunk
            else
                error = abs(result - output)
                fitness -= (1 / (error + 1)) * chunk
            end
            
        catch
            continue
        end
    end

    #@info @sprintf("expr_str: %s\nFitness: %0.2f\n\n", expr_str, fitness)
    indiv.fitness = fitness
end

end

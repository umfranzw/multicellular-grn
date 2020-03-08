module FitnessMod

using IndividualMod
using CellTreeMod
using SymMod
using Printf

function eval(indiv::Individual)
    accuracy = get_accuracy(indiv)
    evolvability = get_evolvability(indiv)
    indiv.fitness = indiv.run.accuracy_weight * accuracy + indiv.run.evolvability_weight * evolvability
end

function get_evolvability(indiv::Individual)
    #number of leaves containing non-terminal symbols
    is_non_term(cell) = cell.sym != nothing && cell.sym.type == SymMod.FcnCall
    num_non_term = 0
    CellTreeMod.traverse(cell -> num_non_term += Int64(is_non_term(cell)), indiv.cell_tree)

    #number of non-app-contributing genes
    app_contrib_genes = ChainGraphMod.get_app_contributing_genes(indiv.chain_graph)
    num_non_contrib = length(indiv.genes) - length(app_contrib_genes)
end

function get_accuracy(inidv::Individual)
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

    fitness
end

end

module FitnessMod

using IndividualMod
using CellTreeMod
using CellMod
using SymMod
using Printf

function eval(indiv::Individual, ea_step::Int64)
    run = indiv.config.run
    completeness = ea_step / run.ea_steps
    
    accuracy = get_accuracy(indiv)
    evolvability = get_evolvability(indiv)
    acc_weight = min(run.initial_acc_weight + completeness * run.weight_shift, run.max_acc_weight)
    ev_weight = max(run.initial_ev_weight - completeness * run.weight_shift, run.min_ev_weight)
    
    indiv.fitness = accuracy * acc_weight + evolvability * ev_weight
end

function get_evolvability(indiv::Individual)
    #number of leaves containing non-terminal symbols
    non_term_fitness = 1.0
    num_non_term = 0
    CellTreeMod.traverse(cell -> num_non_term += Int64(is_non_term(cell)), indiv.cell_tree)
    size = CellTreeMod.size(indiv.cell_tree)
    non_term_fitness -= num_non_term / size

    #number of producing genes
    prod_fitness = 1.0
    num_producing_genes = count(s -> s > 0, indiv.gene_scores)
    prod_fitness -= num_producing_genes / length(indiv.genes)

    (non_term_fitness + prod_fitness) / 2
end

function is_non_term(cell::Cell)
    cell.sym != nothing && cell.sym.type == SymMod.FcnCall
end

function get_accuracy(indiv::Individual)
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
                fitness -= (1 / (1 + error)) * chunk
            end
            
        catch
            continue
        end
    end

    fitness
end

end

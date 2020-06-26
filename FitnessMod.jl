module FitnessMod

using IndividualMod
using CellTreeMod
using CellMod
using SymMod
using Printf

function eval(indiv::Individual, ea_step::Int64)
    # run = indiv.config.run
    # completeness = ea_step / run.ea_steps #low -> high
    # inv_completeness = 1.0 - completeness #high -> low

    # size_fitness = get_size_fitness(indiv, 3)
    # contains_x_fitness = get_contains_x_fitness(indiv)
    # non_term_fitness = get_non_term_fitness(indiv)
    # prod_genes_fitness = get_prod_genes_fitness(indiv)
    # accuracy_fitness = get_accuracy_fitness(indiv)

    # #note:
    # #- (completeness + inv_completeness) = 1
    # #- (completeness * x + inv_completeness * x) = x
    # incr_weight = completeness * 0.6
    # decr_weight = inv_completeness * 0.6
    # stable_weight = 0.4
    
    # fitness = 0.0
    # fitness += (size_fitness + non_term_fitness) * decr_weight / 2
    # fitness += (contains_x_fitness + prod_genes_fitness) * stable_weight / 2
    # fitness += accuracy_fitness * incr_weight

    # indiv.fitness = fitness

    indiv.fitness = get_accuracy_fitness(indiv)
end

function get_size_fitness(indiv::Individual, target_size::Int64)
    size = CellTreeMod.size(indiv.cell_tree)
    
    1.0 - min(size / target_size, 1.0)
end

function get_contains_x_fitness(indiv::Individual)
    contains_x = CellTreeMod.contains_sym(indiv.cell_tree.root, :x)

    Float64(!contains_x)
end

function get_non_term_fitness(indiv::Individual)
    #number of leaves containing non-terminal symbols
    non_term_fitness = 1.0
    num_non_term = 0
    CellTreeMod.traverse(cell -> num_non_term += Int64(is_non_term(cell)), indiv.cell_tree)
    size = CellTreeMod.size(indiv.cell_tree)
    non_term_fitness -= num_non_term / size

    non_term_fitness
end

function get_prod_genes_fitness(indiv::Individual)
    #number of producing genes
    prod_fitness = 1.0
    num_producing_genes = count(s -> s > 0, indiv.gene_scores)
    prod_fitness -= num_producing_genes / length(indiv.genes)

    prod_fitness
end

# function get_evolvability(indiv::Individual)
#     #number of leaves containing non-terminal symbols
#     non_term_fitness = 1.0
#     num_non_term = 0
#     CellTreeMod.traverse(cell -> num_non_term += Int64(is_non_term(cell)), indiv.cell_tree)
#     size = CellTreeMod.size(indiv.cell_tree)
#     non_term_fitness -= num_non_term / size

#     #number of producing genes
#     prod_fitness = 1.0
#     num_producing_genes = count(s -> s > 0, indiv.gene_scores)
#     prod_fitness -= num_producing_genes / length(indiv.genes)

#     (non_term_fitness + prod_fitness) / 2
# end

function is_non_term(cell::Cell)
    cell.sym != nothing && cell.sym.type == SymMod.FcnCall
end

function get_accuracy_fitness(indiv::Individual)
    expr_str = "f(x) = "
    expr_str *= CellTreeMod.to_expr_str(indiv.cell_tree)

    #f(x) = x + 1
    test_data = [(1, 2), (4, 5), (8, 9)]

    #f(x) = x * 2 + 1
    #test_data = [(1, 3), (2, 5), (3, 7)]
    
    fitness = 1.0
    num_data_tests = length(test_data)
    num_extra_tests = 1
    total_tests = num_data_tests + num_extra_tests
    chunk = 1 / total_tests
    num_exceptions = 0

    #input / output tests
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
            num_exceptions += 1
            continue
        end
    end

    #other tests
    #reward a lack of exceptions
    fitness -= (1 / (1 + num_exceptions)) * chunk

    fitness
end

end

module FitnessMod

using IndividualMod
using CellTreeMod
using CellMod
using SymMod
using RegSimInfoMod
using FitnessInfoMod
using Printf
import SettingsMod

function eval_intermediate(indiv::Individual, reg_step::Int64)
    #ideas:
    #reward cycles
    #reward inter-cell signalling
    #reward bind coverage
    #reward prod coverage
    #reward non-terminal leaves (that haven't got a symbol yet)
    #reward the appearance of x
    #reward the duration of production without stagnation
end

#note: there's nothing restricting fitnesses' upper (worst) limit to 1.0!
function eval_final(indiv::Individual, ea_step::Int64)
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
    # fitness += (size_fitness + non_term_fit\ness) * decr_weight / 2
    # fitness += (contains_x_fitness + prod_genes_fitness) * stable_weight / 2
    # fitness += accuracy_fitness * incr_weight

    # indiv.fitness = fitness
    
    indiv.fitness_info = FitnessInfo(
        get_contains_x_fitness(indiv),
        get_contains_fncall_fitness(indiv),
        get_bind_coverage_fitness(indiv),
        get_prod_coverage_fitness(indiv),
        get_divided_fitness(indiv),
        get_altered_sym_prob_fitness(indiv),
        get_genome_len_fitness(indiv),
        get_accuracy_fitness(indiv),
        get_lifetime_fitness(indiv)
    )

    indiv.fitness =
        0.4 +
        0.1 * indiv.fitness_info.contains_fncall +
        0.1 * indiv.fitness_info.divided +
        0.1 * indiv.fitness_info.genome_len +
        0.1 * indiv.fitness_info.prod_coverage + 
        0.1 * indiv.fitness_info.accuracy +
        0.1 * indiv.fitness_info.lifetime
end

function get_lifetime_fitness(indiv::Individual)
    1 - indiv.reg_sim_info.reg_step_count / indiv.config.run.reg_steps
end

function get_genome_len_fitness(indiv::Individual)
    Float64(!(length(indiv.genes) > indiv.config.run.num_initial_genes))
end

function get_altered_sym_prob_fitness(indiv::Individual)
    Float64(!RegSimInfoMod.altered_sym_probs(indiv.reg_sim_info))
end

function get_contains_x_fitness(indiv::Individual)
    get_contains_syms_fitness(indiv, SettingsMod.terms[1:1])
end

function get_contains_fncall_fitness(indiv::Individual)
    get_contains_syms_fitness(indiv, SettingsMod.fcns)
end

function get_contains_syms_fitness(indiv::Individual, syms::Array{Sym, 1})
    contained = CellTreeMod.contains_syms(indiv.cell_tree.root, syms)

    Float64(!contained)
end

function get_bind_coverage_fitness(indiv::Individual)
    bind_coverage = RegSimInfoMod.get_bind_coverage(indiv.reg_sim_info)
    
    1.0 - bind_coverage
end

function get_prod_coverage_fitness(indiv::Individual)
    prod_coverage = RegSimInfoMod.get_prod_coverage(indiv.reg_sim_info)
    
    1.0 - prod_coverage
end

function get_divided_fitness(indiv::Individual)
    Float64(!RegSimInfoMod.divided(indiv.reg_sim_info))
end

function get_accuracy_fitness(indiv::Individual)
    expr_str = "f(x) = "
    expr_str *= CellTreeMod.to_expr_str(indiv.cell_tree)

    #f(x) = x + 1
    test_data = [(1, 2), (4, 5), (8, 9)]

    #f(x) = x * 2 + 1
    #test_data = [(1, 3), (2, 5), (3, 7)]
    
    fitness = 0
    num_data_tests = length(test_data)
    num_extra_tests = 0
    total_tests = num_data_tests + num_extra_tests
    chunk = 1 / total_tests
    #num_exceptions = 0

    #input / output tests
    for (input, output) in test_data
        try
            test_str = "$(expr_str); f($(input))"
            test_expr = Meta.parse(test_str)
            result = Meta.eval(test_expr)

            if result != output
                error = abs(result - output)
                #Notes:
                # 1 / (1 + error) shrinks as error grows (when error is 0, result is 1. When error > 0, result is < 1, range is [0, 1])
                # 1 - 1 / (1 + error) grows as error grows, range is [0, 1]
                # (1 - 1 / (1 + error) * chunk grows as error grows, range is [0, chunk]
                fitness += (1 - 1 / (1 + error)) * chunk
            end
            
        catch
            #num_exceptions += 1
            #consider an exception to be very bad - penalize to the max extent (set error = 1)
            #error = 1
            #fitness += (1 - 1 / (1 + error)) * chunk
            fitness += chunk
            continue
        end
    end

    fitness
end

#------------

function get_size_fitness(indiv::Individual, target_size::Int64)
    size = CellTreeMod.size(indiv.cell_tree)
    
    1.0 - min(size / target_size, 1.0)
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
    num_producing_genes = count(s -> s > 0, indiv.reg_sim_info.produce_count)
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
#     num_producing_genes = count(s -> s > 0, indiv.produce_count)
#     prod_fitness -= num_producing_genes / length(indiv.genes)

#     (non_term_fitness + prod_fitness) / 2
# end

function is_non_term(cell::Cell)
    cell.sym != nothing && cell.sym.type == SymMod.FcnCall
end

end

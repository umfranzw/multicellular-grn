module ProteinAppStepMod

using IndividualMod
using CellTreeMod
using CellMod
using ProteinStoreMod
using ProteinPropsMod
using ProteinAppActionsMod
using ProteinMod
using GeneMod

function run_protein_app(indiv::Individual)
    #we want to visit cells in breadth-first order
    #build an array of the cells (in bfs order) so that as the tree is modified,
    #we don't get messed up by any modifications (eg. new nodes that get added)
    bfs_list = Array{Cell, 1}()
    CellTreeMod.traverse_bf(c -> push!(bfs_list, c), indiv.cell_tree)

    for cell in bfs_list
        run_protein_app_for_cell(indiv.cell_tree, cell, indiv.genes)
    end
end

function run_protein_app_for_cell(tree::CellTree, cell::Cell, genes::Array{Gene, 1})
    #get all proteins (from this cell) that are eligible for application
    app_proteins = ProteinStoreMod.get_by_type(cell.proteins, ProteinPropsMod.Application)

    #build a list of tuples of the form (protein, amount_above_threshold), where each protein has a conc >= protein_app_threshold
    pairs = Array{Tuple{Protein, Float64, Function}, 1}()
    for protein in app_proteins
        #get the appropriate threshold, based on the Protein's Action
        if protein.props.action == ProteinPropsMod.SymProb
            threshold = cell.config.run.sym_prob_threshold
            action_fcn = ProteinAppActionsMod.alter_sym_prob
        elseif protein.props.action == ProteinPropsMod.Divide
            threshold = cell.config.run.cell_division_threshold
            action_fcn = ProteinAppActionsMod.divide
        # elseif protein.props.action == ProteinPropsMod.Sensor
        #     threshold = cell.config.run.sensor_reinforcement_threshold
        #     action_fcn = ProteinAppActionsMod.alter_sensor
        end

        max_conc = maximum(protein.concs)
        excess = max_conc - threshold
        if excess > 0
            push!(pairs, (protein, excess, action_fcn))
        end
    end
    #sort in descending order by sum - we'll apply them in this order
    sort!(pairs; by=p -> p[2], rev=true)

    #apply the proteins
    for pair in pairs
        protein = pair[1]
        action_fcn = pair[3]
        args = AppArgs(tree, cell, genes, protein)

        # protein_str = ProteinPropsMod.to_str(protein.props)
        # println("Applying protein: $(protein_str)")

        action_fcn(args)
    end
end

end

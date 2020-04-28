module ProteinAppActionsMod

using ProteinPropsMod
using ProteinMod
using CellMod
using CellTreeMod
using SymMod
using GeneMod
using ProteinStoreMod
using SymProbsMod
using Printf

export AppArgs

struct AppArgs
    tree::CellTree
    cell::Cell
    genes::Array{Gene, 1}
    app_protein::Protein
end

#note: the methods in this module assume that the app_protein has exceeded the appropriate threshold

function divide(args::AppArgs)
    src_cell = args.cell
    if (length(src_cell.children) < src_cell.config.run.max_children &&
        src_cell.age < src_cell.config.run.division_age_limit &&
        CellTreeMod.size(args.tree) < args.cell.config.run.max_tree_size) #note: save the most expensive check for last
        
        #println("Adding new cell")
        
        max_children = src_cell.config.run.max_children
        num_concs = length(args.genes)
        chunk_index = args.app_protein.props.arg % max_children #in [0, max_children - 1]
        truncated_chunk_size = num_concs ÷ max_children
        chunk_start = chunk_index * truncated_chunk_size + 1
        #last chunk gets the extra concs (if there are any)
        chunk_end = (chunk_index == max_children - 1) ? num_concs : chunk_start + truncated_chunk_size
        chunk_range = chunk_start:chunk_end

        new_cell = Cell(src_cell.config, args.genes, -1) #age is -1 since it'll be incremented at end of this reg step, and we want this cell to be younger than its parent even if the division occurs on reg step 1
        CellMod.add_parent(new_cell, src_cell)
        

        for protein in ProteinStoreMod.get_all(src_cell.proteins)
            new_protein = Protein(src_cell.config, deepcopy(protein.props), false, false, num_concs, pointer_from_objref(src_cell))
            #child gains half of the (magnitude of the) parent's concs in the chunk range
            protein.concs[chunk_range] = protein.concs[chunk_range] ./ 2
            new_protein.concs[chunk_range] = protein.concs[chunk_range]

            ProteinStoreMod.insert(new_cell.proteins, new_protein)
        end
    end

    #remove the app_protein, since it's done its job
    ProteinStoreMod.remove(src_cell.proteins, args.app_protein)
end

function alter_sym_prob(args::AppArgs)
    cell = args.cell
    if cell.sym == nothing #if symbol has not yet been fixed
        #println("Altering Sym Prob")
        age_factor = 1.0 - cell.age / cell.config.run.reg_steps
        
        sym_index = Int64(args.app_protein.props.arg % length(cell.probs)) + 1
        sign = args.app_protein.props.fcn == ProteinPropsMod.Inhibit ? -1 : 1
        
        max_excess = 1.0 - cell.config.run.sym_prob_threshold #max possible excess
        scale_factor = 1.0 / max_excess
        excess = maximum(args.app_protein.concs) - cell.config.run.sym_prob_threshold #should be positive, given that this method has been called
        delta = sign * excess * scale_factor * age_factor
        SymProbsMod.alter_prob(cell.probs, sym_index, delta)
    end

    #remove the app_protein, since it's done its job
    ProteinStoreMod.remove(cell.proteins, args.app_protein)
end

function alter_sensor(args::AppArgs)
    #println("Altering sensor")
    
    cell = args.cell
    num_concs = length(args.app_protein.concs)
    loc = args.app_protein.props.loc

    #sign = args.app_protein.props.fcn == ProteinPropsMod.Inhibit ? -1 : 1
    max_excess = 1.0 - cell.config.run.sensor_reinforcement_threshold #max possible excess
    scale_factor = 1.0 / max_excess
    excess = max.(args.app_protein.concs .- cell.config.run.sensor_reinforcement_threshold, 0.0)
    #delta = sign * excess .* scale_factor
    delta = excess .* scale_factor
    CellMod.adjust_sensor(cell, loc, delta)

    ProteinStoreMod.remove(cell.proteins, args.app_protein)
end

end

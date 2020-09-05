module ProteinAppActionsMod

using ProteinPropsMod
using ProteinMod
using CellMod
using CellTreeMod
using SymMod
using GeneMod
using ProteinStoreMod
using SymProbsMod
using RegSimInfoMod
using Printf

export AppArgs

struct AppArgs
    info::RegSimInfo
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
        new_cell = Cell(src_cell.config, args.genes, -1) #age is -1 since it'll be incremented at end of this reg step, and we want this cell to be younger than its parent even if the division occurs on reg step 1
        CellMod.add_parent(new_cell, src_cell)

        # num_concs = length(args.genes)
        # for protein in ProteinStoreMod.get_by_type(src_cell.proteins, ProteinPropsMod.Internal)
        #     new_protein = Protein(src_cell.config, deepcopy(protein.props), false, false, num_concs, src_cell.id)
        #     #child gains half of the (magnitude of the) parent's concs
        #     new_protein.concs = protein.concs ./ 2

        #     ProteinStoreMod.insert(new_cell.proteins, new_protein)
        # end
        args.info.division_count += 1
    end
end

function alter_sym_prob(args::AppArgs)
    cell = args.cell
    if cell.sym == nothing #if symbol has not yet been fixed
        # println("Altering Sym Prob")
        # println("cell ID: $(cell.id)")
        # show(args.app_protein)
        
        age_factor = 1.0 - cell.age / cell.config.run.reg_steps
        
        sym_index = abs(Int64(args.app_protein.props.arg)) % length(cell.probs) + 1
            
        # if sym_index < 0
        #     println(args.app_protein.props.arg)
        #     println(abs(args.app_protein.props.arg))
        #     println(abs(args.app_protein.props.arg) % length(cell.probs))
        #     println(Int64(abs(args.app_protein.props.arg) % length(cell.probs)))
        #     println(Int64(abs(args.app_protein.props.arg) % length(cell.probs)) + 1)
        # end
            
        sign = Int64(ProteinPropsMod.get_fcn(args.app_protein.props))
        
        max_excess = 1.0 - cell.config.run.sym_prob_threshold #max possible excess
        scale_factor = 1.0 / max_excess
        excess = maximum(args.app_protein.concs) - cell.config.run.sym_prob_threshold #should be positive, given that this method has been called
        delta = sign * excess * scale_factor * age_factor
        SymProbsMod.alter_prob(cell.probs, sym_index, delta)
        args.info.alter_sym_prob_count += 1
    end
end

function alter_sensor(args::AppArgs)
    #println("Altering sensor")
    
    cell = args.cell
    num_concs = length(args.app_protein.concs)
    max_locs = 3 + cell.config.run.max_children
    loc = Int8(abs(args.app_protein.props.arg) % max_locs)

    max_excess = 1.0 - cell.config.run.sensor_reinforcement_threshold #max possible excess
    scale_factor = 1.0 / max_excess
    excess = max.(args.app_protein.concs .- cell.config.run.sensor_reinforcement_threshold, 0.0)
    #delta = sign * excess .* scale_factor
    delta = excess .* scale_factor
    CellMod.adjust_sensor(cell, loc, delta)
end

end

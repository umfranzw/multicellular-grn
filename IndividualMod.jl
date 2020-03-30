module IndividualMod

using RunMod
using GeneMod
using GeneStateMod
using ProteinMod
using ProteinPropsMod
using ProteinAppActionsMod
using CellMod
using CellTreeMod
using ProteinStoreMod
using SymMod
using DiffusionMod
using CellTreeMod
using MiscUtilsMod
using Printf
using ChainGraphMod

import Random
import RandUtilsMod
import Base.show
import Random

export Individual,
    rand_init, run_bind

mutable struct Individual
    config::Config
    genes::Array{Gene, 1}
    cell_tree::CellTree
    initial_cell_proteins::Array{Protein, 1}
    #note: this is a value in [0.0, 1.0], where 0.0 is optimal
    fitness::Float64
    chain_graph::ChainGraph
end

function show(io::IO, indiv::Individual, ilevel::Int64=0)
    iprintln(io, "Individual:", ilevel)
    iprintln(io, "Genes:", ilevel + 1)
    map(g -> GeneMod.show(io, g, ilevel + 2), indiv.genes)

    iprintln(io, "cell_tree:", ilevel + 1)
    iprintln(io, indiv.cell_tree, ilevel + 2)

    iprintln(io, "initial_cell_proteins:", ilevel + 1)
    map(p -> ProteinMod.show(io, p, ilevel + 2), indiv.initial_cell_proteins)

    iprintln(io, "fitness: $(indiv.fitness)", ilevel + 1)
end

function rand_init(run::Run, seed::UInt64)
    rng = Random.MersenneTwister(seed)
    config = Config(run, rng)
    
    genes = map(i -> GeneMod.rand_init(config, i), 1:config.run.num_initial_genes)
    root_cell = Cell(config, genes)
    cell_tree = CellTree(root_cell)
    
    initial_proteins = Array{Protein, 1}()
    
    for i in 1:config.run.num_initial_proteins
        #all initial proteins should have type=Internal
        props = ProteinPropsMod.rand_init(config, type=[ProteinPropsMod.Internal]])

        #note: it is possible that not all initial proteins in the array are unique. That's ok, since they'll be subject to evolution.
        protein = Protein(config, props, true, true, length(root_cell.gene_states), pointer_from_objref(root_cell))
        push!(initial_proteins, protein)
    end
    
    indiv = Individual(config, genes, cell_tree, initial_proteins, 1.0, ChainGraph())
    CellMod.insert_initial_proteins(root_cell, indiv.initial_cell_proteins)

    indiv
end

#resets everything to the way it was before the reg sim (so
#we can run the reg sim again on the next ea_step)
function reset(indiv::Individual)
    #just re-initialize the cell (this discards the rest of the tree, along with any protein bindings)
    indiv.cell_tree.root = Cell(indiv.config, indiv.genes)
    CellMod.insert_initial_proteins(indiv.cell_tree.root, indiv.initial_cell_proteins)
end

function run_protein_app(indiv::Individual)
    #we want to visit cells in breadth-first order
    #build an array of the cells (in bfs order) so that as the tree is modified,
    #we don't get messed up by any modifications (eg. new nodes that get added)
    bfs_list = Array{Cell, 1}()
    CellTreeMod.traverse_bf(c -> push!(bfs_list, c), indiv.cell_tree)

    deleted_cells = Set{Cell}()
    for cell in bfs_list
        if cell.energy > indiv.config.run.cell_energy_threshold && cell ∉ deleted_cells
            deleted = run_protein_app_for_cell(indiv.cell_tree, cell, indiv.genes, indiv.initial_cell_proteins)
            deleted_cells = union(deleted_cells, deleted...)
        end
    end
end

function run_protein_app_for_cell(tree::CellTree, cell::Cell, genes::Array{Gene, 1}, initial_proteins::Array{Protein, 1})
    #get all proteins (from this cell) that are eligible for application
    app_proteins = ProteinStoreMod.get_by_type(cell.proteins, ProteinPropsMod.App)

    #build a list of tuples of the form (protein, amount_above_threshold), where each protein has a conc >= protein_app_threshold
    pairs = Array{Tuple{Protein, Float64}, 1}()
    for protein in app_proteins
        comp = repeat([cell.config.run.protein_app_threshold], length(protein.concs))
        result = protein.concs .- comp
        overflow = max(result...)
        if overflow >= 0
            push!(pairs, (protein, overflow))
        end
        
        # conc_sum = sum(protein.concs)
        # if conc_sum >= cell.config.run.protein_app_threshold
        #     push!(pairs, (protein, conc_sum))
        # end
    end
    #sort in descending order by sum - we'll apply they in this order
    sort!(pairs; by=p -> p[2], rev=true)

    #apply the proteins
    deleted_cells = Set{Cell}()
    for pair in pairs
        protein = pair[1]
        deleted = ProteinAppActionsMod.run_app_action(tree, cell, genes, initial_proteins, protein)
        deleted_cells = union(deleted_cells, deleted...)
    end

    deleted_cells
end

function update_chains(indiv::Individual, reg_step::Int64)
    ChainGraphMod.append_for_tree(indiv.chain_graph, indiv.cell_tree, reg_step)
end

function run_bind(indiv::Individual)
    CellTreeMod.traverse(cell -> run_bind_for_cell(indiv, cell), indiv.cell_tree)
end

function run_produce(indiv::Individual)
    CellTreeMod.traverse(cell -> run_produce_for_cell(indiv, cell), indiv.cell_tree)
end

function run_diffuse(indiv::Individual)
    DiffusionMod.diffuse_intra_cell_proteins(indiv.cell_tree)
    DiffusionMod.diffuse_inter_cell_proteins(indiv.cell_tree)
end

function run_decay(indiv::Individual)
    CellTreeMod.traverse(cell -> run_decay_for_cell(cell), indiv.cell_tree)
end

function run_neighbour_comm(indiv::Individual)
    CellTreeMod.traverse(cell -> run_neighbour_comm_for_cell(cell), indiv.cell_tree)
end

function run_neighbour_comm_for_cell(cell::Cell)
    for src_loc in instances(ProteinProps.ProteinLoc)
        neighbour = get_neighbour(src_loc)
        sensor_amount = cell.sensors[src_loc]

        #get neighbour's neighbour proteins for opposite loc (the ones that are being sent to this cell)
        opposite_loc = get_opposite_loc(src_loc)
        neighbour_proteins = ProteinStoreMod.get_neighbour_proteins_by_loc(neighbour.proteins, opposite_loc, ptr_from_objectref(neighbour))

        #compute the max amount of each protein that we can accept
        accept_amount = sensor_amount / length(neighbour_proteins)
        for neighbour_protein in neighbour_proteins
            transfer_amount = min.(protein.concs, accept_amount)
            cell.sensors[src_loc] -= transfer_amount
            neighbour_protein.concs -= amount
            
            #add the neighbour protein into the source cell, flipping the loc
            dest_props = ProteinProps(
                cell.config,
                neighbour_protein.type,
                neighbour_protein.fcn,
                neighbour_protein.action,
                src_loc, #flip the loc
                neighbour_protein.arg
            )

            dest_protein = ProteinStoreMod.get(cell.proteins, props)
            if dest_protein == nothing
                dest_protein = Protein(cell.config, dest_props, false, false, length(cell.genes), ptr_from_objectref(neighbour))
            end
            dest_protein.concs = clamp.(dest_protein.concs + transfer_amount, 0.0, 1.0)

            #note: if the neighbour_protein was completely tranfered, the decay step will delete it later
        end
    end
end

function run_decay_for_cell(cell::Cell)
    proteins = ProteinStoreMod.get_all(cell.proteins)
    for protein in proteins
        #decrease the concentration using deca_rate
        protein.concs = 
            protein.concs = max.(protein.concs .- protein.concs * cell.config.run.decay_rate, zeros(length(protein.concs)))

        #remove any proteins that have decayed below the allowable threshold
        if all(c -> c < cell.config.run.protein_deletion_threshold, protein.concs)
            #clear any bindings that this protein has
            for gs in cell.gene_states
                for i in 1:length(gs.bindings)
                    if gs.bindings[i] == protein
                        gs.bindings[i] = nothing
                    end
                end
            end

            #remove it from the cell
            ProteinStoreMod.remove(cell.proteins, protein)
        end
    end
end

function run_produce_for_cell(indiv::Individual, cell::Cell)
    for gene_index in 1:length(cell.gene_states)
        gene_state = cell.gene_states[gene_index]
        gene = indiv.genes[gene_index]
        #note: rates is an array of the form [(prod_site_index, rate), ...]
        rates = GeneStateMod.get_prod_rates(gene_state)

        for (prod_index, rate) in rates
            run_produce_for_site(cell, gene_index, prod_index, rate)
        end
    end
end

function run_produce_for_site(cell::Cell, gene_index::Int64, prod_index::Int64, rate::Float64)
    #get the props for the protein that will be produced
    gene = cell.gene_states[gene_index].gene
    props = gene.prod_sites[prod_index]
    #check if protein already exists in this cell's store
    protein = ProteinStoreMod.get(cell.proteins, props)
    #if not, create and insert it
    if protein == nothing
        #note: protein will be initialized with conc values of zero
        #@info @sprintf("Produced protein: %s", props)
        protein = Protein(cell.config, props, false, false, length(cell.gene_states), pointer_from_objref(cell))
        ProteinStoreMod.insert(cell.proteins, protein)
    end
    
    #increment the conc using the rate
    #note: this will only increment the conc directly over the gene
    #the diffusion will spread this out later
    protein.concs[gene_index] = min(protein.concs[gene_index] + rate, 1.0)
end

function run_bind_for_cell(indiv::Individual, cell::Cell)
    for gene_index in 1:length(cell.gene_states)
        gene = indiv.genes[gene_index]

        for i in 1:length(gene.bind_sites)
            eligible_proteins = get_bind_eligible_proteins_for_site(cell, gene, bind_site_index)
            run_bind_for_site(cell.gene_states[gene_index], gene_index, site_index, elibible_proteins)
        end
    end
end

function get_bind_eligible_proteins_for_site(cell::Cell, gene::Gene, site_index::Int64)
    #Binding logic:
    #  -protein's type must match site's type
    #  -protein's conc must be >= site's threshold
    #  -if the type is neighbour or diffusion, the protein's src_cell_ptr must not be equal to this cell ptr (to prevent self-binding)
    #  -protein's action must match site's action
    #  -protein's loc must match site's loc

    site = gene.bind_sites[site_index]
    eligible_proteins = Array{Protein, 1}()
    for protein in values(ProteinStoreMod.get_by_type(cell.proteins, site.type))
        eligible = protein.concs[gene_index] >= site.bind_threshold
        if site.type == ProteinPropsMod.Neighbour || site.type == ProteinPropsMod.Diffusion
            eligible = eligible && protein.src_cell_ptr != pointer_from_objref(cell)
        end
        
        eligible = eligible && protein.props.action == site.action
        eligible = eligible && protein.props.loc == site.loc
        
        if eligible
            push!(eligible_proteins, protein)
        end
    end

    eligible_proteins
end

function run_bind_for_site(gs::GeneState, gene_index::Int64, site_index::Int64, eligible_proteins::Array{Protein, 1})
    if length(eligible_proteins) > 0
        #use roulette wheel style selection to pick the protein
        conc_sum = foldl((s, p) -> s + p.concs[gene_index], eligible_proteins; init=0.0)
        next_bound = 0.0
        wheel = []
        for protein in eligible_proteins
            next_bound += protein.concs[gene_index] / conc_sum
            push!(wheel, next_bound)
        end

        wheel[end] = 1.0 #just in case we've got some floating point error

        sel_index = 1
        r = RandUtilsMod.rand_float(gs.config) #random value in [0, 1)
        while r >= wheel[sel_index]
            sel_index += 1
        end

        sel_protein = eligible_proteins[sel_index]

        #@info @sprintf("%s binding to site %s", sel_protein, site_type)
        GeneStateMod.bind(gs, sel_protein, site_index)

    else
        state = GeneStateMod.get_binding(gs, site_index)
        if state != nothing
            #@info @sprintf("Protein unbinding")
            GeneStateMod.unbind(gs, site_index)
        end
    end
end

end

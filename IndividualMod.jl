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
using SymProbsMod
using Printf

import Random
import RandUtilsMod
import Base.show
import Random

export Individual,
    rand_init, run_bind

initial_threshold = 0.1
initial_consum_rate = 0.01

mutable struct Individual
    config::Config
    genes::Array{Gene, 1}
    cell_tree::CellTree
    initial_cell_proteins::Array{Protein, 1}
    #note: this is a value in [0.0, 1.0], where 0.0 is optimal
    fitness::Float64
    gene_scores::Array{Int64, 1}
end

function rand_init(run::Run, seed::UInt64)
    rng = Random.MersenneTwister(seed)
    config = Config(run, rng)

    genes, initial_protein_props = make_initial_genes(config)
    #genes = map(i -> GeneMod.rand_init(config, i, [ProteinPropsMod.Internal], [GeneMod.Id]), 1:config.run.num_initial_genes)
    
    root_cell = Cell(config, genes)
    cell_tree = CellTree(root_cell)
    
    initial_proteins = make_initial_proteins(config, initial_protein_props, length(genes), root_cell)
    
    indiv = Individual(config, genes, cell_tree, initial_proteins, 1.0, zeros(Int64, length(genes)))
    CellMod.insert_initial_proteins(root_cell, indiv.initial_cell_proteins)

    indiv
end

# function make_initial_genes(config::Config)
#     genes = Array{Gene, 1}()
#     initial_protein_props = Array{ProteinProps, 1}()
#     genome_index = 0

#     genome_index += 1
#     g1_bind_site = GeneMod.rand_bind_site(
#         config,
#         type=[ProteinPropsMod.Internal],
#         threshold=[IndividualMod.initial_threshold],
#         consum_rate=[IndividualMod.initial_consum_rate]
#     )

#     g1_prod_site = GeneMod.rand_prod_site(
#         config,
#         type=[ProteinPropsMod.Internal],
#         fcn=[ProteinPropsMod.Activate]
#     )
#     g1 = pop_remaining_sites(config, genome_index, g1_bind_site, g1_prod_site)
#     push!(genes, g1)

#     props = ProteinPropsMod.rand_init(
#         config,
#         type=[g1_bind_site.type],
#         fcn=[ProteinPropsMod.Activate],
#         action=[g1_bind_site.action],
#         loc=[g1_bind_site.loc]
#     )
#     push!(initial_protein_props, props)

#     (genes, initial_protein_props)
# end

function make_initial_genes(config::Config)
    genes = Array{Gene, 1}()
    initial_protein_props = Array{ProteinProps, 1}()
    genome_index = -1

    #ab
    ab_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ProteinPropsMod.Internal],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    ab_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Internal],
        fcn=[ProteinPropsMod.Activate]
    )
    ab = pop_remaining_sites(config, genome_index, ab_bind_site, ab_prod_site)

    props = ProteinPropsMod.rand_init(
        config,
        type=[ab_bind_site.type],
        fcn=[ProteinPropsMod.Activate],
        action=[ab_bind_site.action],
        loc=[ab_bind_site.loc]
    )
    push!(initial_protein_props, props)

    #bc
    bc_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ab_prod_site.type],
        action=[ab_prod_site.action],
        loc=[ab_prod_site.loc],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    bc_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Application],
        action=[ProteinPropsMod.Divide]
    )
    bc = pop_remaining_sites(config, genome_index, bc_bind_site, bc_prod_site)

    #ba
    ba_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ab_prod_site.type],
        action=[ab_prod_site.action],
        loc=[ab_prod_site.loc],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    ba_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ab_bind_site.type],
        fcn=[ProteinPropsMod.Activate],
        action=[ab_bind_site.action],
        loc=[ab_bind_site.loc]
    )
    ba = pop_remaining_sites(config, genome_index, ba_bind_site, ba_prod_site)

    #de
    de_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ProteinPropsMod.Internal],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    plus_index = findfirst(i -> SymProbsMod.index_to_sym[i].val == :+, 1:length(SymProbsMod.index_to_sym)) - 1
    de_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Application],
        fcn=[ProteinPropsMod.Activate],
        action=[ProteinPropsMod.SymProb],
        arg=[UInt8(plus_index)]
    )
    de = pop_remaining_sites(config, genome_index, de_bind_site, de_prod_site)
    
    props = ProteinPropsMod.rand_init(
        config,
        type=[de_bind_site.type],
        fcn=[ProteinPropsMod.Activate],
        action=[de_bind_site.action],
        loc=[de_bind_site.loc]
    )
    push!(initial_protein_props, props)
    
    #bf
    bf_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ab_prod_site.type],
        action=[ab_prod_site.action],
        loc=[ab_prod_site.loc],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    bf_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Neighbour],
        fcn=[ProteinPropsMod.Activate],
        loc=[ProteinPropsMod.BottomRight]
    )
    bf = pop_remaining_sites(config, genome_index, bf_bind_site, bf_prod_site)

    #bg
    bg_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ab_prod_site.type],
        action=[ab_prod_site.action],
        loc=[ab_prod_site.loc],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    bg_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Neighbour],
        fcn=[ProteinPropsMod.Activate],
        loc=[ProteinPropsMod.BottomLeft]
    )
    bg = pop_remaining_sites(config, genome_index, bg_bind_site, bg_prod_site)

    #gh
    gh_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ProteinPropsMod.Internal],
        action=[bg_prod_site.action],
        loc=[bg_prod_site.loc],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    x_index = findfirst(i -> SymProbsMod.index_to_sym[i].val == :x, 1:length(SymProbsMod.index_to_sym)) - 1
    gh_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Application],
        fcn=[ProteinPropsMod.Activate],
        action=[ProteinPropsMod.SymProb],
        arg=[UInt8(x_index)]
    )
    gh = pop_remaining_sites(config, genome_index, gh_bind_site, gh_prod_site)

    #fnb
    fnb_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ProteinPropsMod.Internal],
        action=[bf_prod_site.action],
        loc=[bf_prod_site.loc],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    fnb_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Neighbour],
        fcn=[ProteinPropsMod.Inhibit],
        loc=[ProteinPropsMod.Top]
    )
    fnb = pop_remaining_sites(config, genome_index, fnb_bind_site, fnb_prod_site)

    #fi
    fi_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ProteinPropsMod.Internal],
        action=[bf_prod_site.action],
        loc=[bf_prod_site.loc],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    one_index = findfirst(i -> SymProbsMod.index_to_sym[i].val == :1, 1:length(SymProbsMod.index_to_sym)) - 1
    fi_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Application],
        fcn=[ProteinPropsMod.Activate],
        action=[ProteinPropsMod.SymProb],
        arg=[UInt8(one_index)]
    )
fi = pop_remaining_sites(config, genome_index, fi_bind_site, fi_prod_site)
push!(genes, gh, bg, ba, ab, bc, fnb, bf, fi, de)
for i in 1:length(genes)
    genes[i].genome_index = i
end

    (genes, initial_protein_props)
end

function pop_remaining_sites(config::Config, genome_index::Int64, first_bind_site::BindSite, first_prod_site::ProdSite)
    bind_sites = Array{BindSite, 1}()
    prod_sites = Array{ProdSite, 1}()
    push!(bind_sites, first_bind_site)
    push!(prod_sites, first_prod_site)
    
    #randomly init the remaining bind and prod sites
    for j in 2:config.run.bind_sites_per_gene
        push!(bind_sites, GeneMod.rand_bind_site(config))
        push!(prod_sites, GeneMod.rand_prod_site(config))
    end

    Gene(config, genome_index, GeneMod.Id, bind_sites, prod_sites)
end

function make_initial_proteins(config::Config, initial_props::Array{ProteinProps, 1}, num_concs::Int64, root_cell::Cell)
    proteins = Array{Protein, 1}()
    half = (1.0 - IndividualMod.initial_threshold) / 2
    for props in initial_props
        #note: it is possible that not all initial proteins in the array are unique. That's ok, since they'll be subject to evolution.
        protein = Protein(config, props, false, true, length(root_cell.gene_states), root_cell.id)
        protein.concs = RandUtilsMod.rand_floats(config, IndividualMod.initial_threshold + half, 1.0, num_concs)
        push!(proteins, protein)
    end

    proteins
end

#resets everything to the way it was before the reg sim (so
#we can run the reg sim again on the next ea_step)
function reset_cell_tree(indiv::Individual)
    #just re-initialize the cell (this discards the rest of the tree, along with any protein bindings)
    indiv.cell_tree.root = Cell(indiv.config, indiv.genes)
    CellMod.insert_initial_proteins(indiv.cell_tree.root, indiv.initial_cell_proteins)
end

function reset_gene_scores(indiv::Individual)
    indiv.gene_scores .= 0
end

function extend_sensors(indiv::Individual, index::Int64)
    CellTreeMod.traverse(cell -> CellMod.extend_sensors(cell, index), indiv.cell_tree)
end

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

function show(io::IO, indiv::Individual, ilevel::Int64=0)
    iprintln(io, "Individual:", ilevel)
    iprintln(io, "genes:", ilevel + 1)
    foreach(g -> GeneMod.show(io, g, ilevel + 2), indiv.genes)

    iprintln(io, "gene_scores:", ilevel + 1)
    iprintln(io, join(indiv.gene_scores, ", "), ilevel + 2)
             
    iprintln(io, "cell_tree:", ilevel + 1)
    iprintln(io, indiv.cell_tree, ilevel + 2)

    iprintln(io, "initial_cell_proteins:", ilevel + 1)
    foreach(p -> ProteinMod.show(io, p, ilevel + 2), indiv.initial_cell_proteins)

    iprintln(io, "fitness: $(indiv.fitness)", ilevel + 1)
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
        elseif protein.props.action == ProteinPropsMod.Sensor
            threshold = cell.config.run.sensor_reinforcement_threshold
            action_fcn = ProteinAppActionsMod.alter_sensor
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

function run_bind(indiv::Individual)
    CellTreeMod.traverse(cell -> run_bind_for_cell(indiv, cell), indiv.cell_tree)
end

function run_produce(indiv::Individual)
    CellTreeMod.traverse(cell -> run_produce_for_cell(indiv, cell), indiv.cell_tree)
end

function run_binding_consum(indiv::Individual)
    CellTreeMod.traverse(run_binding_consum_for_cell, indiv.cell_tree)
end

function run_diffuse(indiv::Individual)
    DiffusionMod.diffuse_intra_cell_proteins(indiv.cell_tree)
    DiffusionMod.diffuse_inter_cell_proteins(indiv.cell_tree)
end

function run_decay(indiv::Individual)
    CellTreeMod.traverse(run_decay_for_cell, indiv.cell_tree)
end

function run_age(indiv::Individual)
    CellTreeMod.traverse(cell -> cell.age += 1, indiv.cell_tree)
end

function run_fix_syms(indiv::Individual)
    CellTreeMod.traverse(CellMod.fix_sym, indiv.cell_tree)
end

function run_neighbour_comm(indiv::Individual)
    info = TreeInfo(indiv.cell_tree)
    CellTreeMod.traverse(cell -> run_neighbour_comm_for_cell(cell, info), indiv.cell_tree)
end

function get_neighbours_at(cell::Cell, info::TreeInfo, loc::ProteinPropsMod.ProteinLoc)
    neighbours = Array{Cell, 1}()
    
    if loc == ProteinPropsMod.Top
        if cell.parent != nothing
            push!(neighbours, cell.parent)
        end
        
    elseif loc == ProteinPropsMod.BottomLeft
        num_children = length(cell.children)
        if num_children > 0
            range = 1:max(num_children ÷ 2, 1)
            append!(neighbours, cell.children[range])
        end

    elseif loc == ProteinPropsMod.BottomRight
        num_children = length(cell.children)
        if num_children > 1
            range = (num_children ÷ 2 + 1):num_children
            append!(neighbours, cell.children[range])
        end

    elseif loc == ProteinPropsMod.Left || loc == ProteinPropsMod.Right
        level = info.cell_to_level[cell]
        row = info.level_to_cell[level]
        cell_index = findall(c -> c == cell, row)[1]

        if loc == ProteinPropsMod.Left
            if cell_index > 1
                push!(neighbours, row[cell_index - 1])
            end
        else
            if cell_index < length(row)
                push!(neighbours, row[cell_index + 1])
            end
        end
    end

    neighbours
end

function run_neighbour_comm_for_cell(cell::Cell, info::TreeInfo)
    for src_loc in instances(ProteinPropsMod.ProteinLoc)
        neighbours = get_neighbours_at(cell, info, src_loc) #in the case that src_loc == BotttomLeft or src_loc == BottomRight, we may have mulitple neighbours (the children)
        for neighbour in neighbours
            sensor_amount = cell.sensors[src_loc] / length(neighbours) #allocate an equal amount to each neighbour

            #get neighbour's neighbour proteins for opposite locs (the ones that are being sent to this cell)
            #note that there may be more than one opposite loc (eg. if src_loc == top, we have bottomleft and bottomright)
            opposite_locs = ProteinPropsMod.get_opposite_locs(src_loc)
            neighbour_proteins = Array{Protein, 1}()
            for opp_loc in opposite_locs
                append!(neighbour_proteins, ProteinStoreMod.get_neighbour_proteins_by_loc(neighbour.proteins, opp_loc, neighbour.id))
            end
            
            #further filter down to the ones that originated in the neighbour cell (it's possible they may have been transferred in from another neighbour)
            neighbour_proteins = filter(p -> p.src_cell_id == neighbour.id, neighbour_proteins)

            #compute the max amount of each protein that we can accept
            accept_amount = sensor_amount / length(neighbour_proteins)
            for neighbour_protein in neighbour_proteins
                transfer_amount = min.(neighbour_protein.concs, accept_amount)
                cell.sensors[src_loc] -= transfer_amount
                neighbour_protein.concs -= transfer_amount
                
                #add the neighbour protein into the source cell
                dest_protein = ProteinStoreMod.get(cell.proteins, neighbour_protein.props)
                if dest_protein == nothing
                    if ProteinStoreMod.num_proteins(cell.proteins) < cell.config.run.max_proteins_per_cell
                        dest_protein = Protein(cell.config, deepcopy(neighbour_protein.props), false, false, length(cell.gene_states), neighbour.id)
                        ProteinStoreMod.insert(cell.proteins, dest_protein)
                    end
                end
                if dest_protein != nothing
                    dest_protein.concs = clamp.(dest_protein.concs + transfer_amount, 0.0, 1.0)
                end

                #note: if the neighbour_protein was completely tranfered, the decay step will delete it later
            end
        end
    end
end

function run_decay_for_cell(cell::Cell)
    proteins = ProteinStoreMod.get_all(cell.proteins)
    for protein in proteins
        #decrease the concentration
        protein.concs = max.(protein.concs .- cell.config.run.decay_rate, zeros(length(protein.concs)))

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

    #run decay on the cell sensors
    # for loc in keys(cell.sensors)
    #     cell.sensors[loc] = max.(cell.sensors[loc] .- cell.config.run.decay_rate, zeros(length(cell.sensors[loc])))
    #     #note: sensor proteins are never removed
    # end
end

function run_produce_for_cell(indiv::Individual, cell::Cell)
    for gene_index in 1:length(cell.gene_states)
        gene_state = cell.gene_states[gene_index]
        gene = indiv.genes[gene_index]
        #note: rates is an array of the form [(prod_site_index, rate), ...]
        rates = GeneStateMod.get_prod_rates(gene_state)

        for (prod_index, rate) in rates
            run_produce_for_site(cell, gene_index, prod_index, rate)
            indiv.gene_scores[gene_index] += 1
        end
    end
end

function run_produce_for_site(cell::Cell, gene_index::Int64, prod_index::Int64, rate::Float64)
    #get the props for the protein that will be produced
    gene = cell.gene_states[gene_index].gene
    prod_site = gene.prod_sites[prod_index]
    props = ProteinProps(
        prod_site.type,
        prod_site.fcn,
        prod_site.action,
        prod_site.loc,
        prod_site.arg
    )
    
    #check if protein already exists in this cell's store
    protein = ProteinStoreMod.get(cell.proteins, props)
    #if not, create and insert it
    if protein == nothing
        if ProteinStoreMod.num_proteins(cell.proteins) < cell.config.run.max_proteins_per_cell
            #note: protein will be initialized with conc values of zero
            #@info @sprintf("Produced protein: %s", props)
            #note: remember to deepcopy the props, so that if this prod site is mutated, the protein's props don't also
            protein = Protein(cell.config, deepcopy(props), false, false, length(cell.gene_states), cell.id)
            ProteinStoreMod.insert(cell.proteins, protein)
        end
    end

    if protein != nothing
        #increment the conc using the rate
        #note: this will only increment the conc directly over the gene
        #the diffusion will spread this out later
        protein.concs[gene_index] = min(protein.concs[gene_index] + rate, 1.0)
    end
end

function run_binding_consum_for_cell(cell::Cell)
    foreach(GeneStateMod.run_binding_consum, cell.gene_states)
end

function run_bind_for_cell(indiv::Individual, cell::Cell)
    for gene_index in 1:length(cell.gene_states)
        gene = indiv.genes[gene_index]

        for site_index in 1:length(gene.bind_sites)
            eligible_proteins = get_bind_eligible_proteins_for_site(cell, gene, site_index)
            run_bind_for_site(indiv.config, cell.gene_states[gene_index], gene_index, site_index, eligible_proteins)
        end
    end
end

function get_bind_eligible_proteins_for_site(cell::Cell, gene::Gene, site_index::Int64)
    #Binding logic:
    # For proteins of type Internal:
    # -protein's conc must be >= site's threshold
    # -protein type must match site type
    # -protein fcn is unrestricted
    # -protein action must match site action
    # -protein loc must match site loc
    # -protein arg is unrestricted

    #For proteins of type Neighbour:
    # -same restrictions as type Internal, except:
    # -bind site type must be internal
    # -protein's src_cell_id must not point to the current cell (to prevent self-binding)

    #For proteins of type Diffusion:
    # -same restrictions as type Neighbour

    #Proteins of type Application do not bind
    
    site = gene.bind_sites[site_index]
    eligible_proteins = Array{Protein, 1}()
    
    #sites with type Internal can accept proteins of types: Internal, Neighbour, or Diffusion
    #note: all bind sites are now of type Internal
    search_proteins = ProteinStoreMod.get_by_types(cell.proteins, Set{ProteinPropsMod.ProteinType}( (ProteinPropsMod.Internal, ProteinPropsMod.Neighbour, ProteinPropsMod.Diffusion) ))
    
    for protein in search_proteins
        eligible = protein.concs[gene.genome_index] >= site.threshold && protein.props.action == site.action && protein.props.loc == site.loc
        if eligible && (protein.props.type == ProteinPropsMod.Neighbour || protein.props.type == ProteinPropsMod.Diffusion)
            eligible = protein.src_cell_id != cell.id
        end
        
        if eligible 
            push!(eligible_proteins, protein)
        end
    end

    eligible_proteins
end

function run_bind_for_site(config::Config, gs::GeneState, gene_index::Int64, site_index::Int64, eligible_proteins::Array{Protein, 1})
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
        r = RandUtilsMod.rand_float(config) #random value in [0, 1)
        while r >= wheel[sel_index]
            sel_index += 1
        end

        sel_protein = eligible_proteins[sel_index]

        #@info @sprintf("%s binding to site %s", sel_protein, site_index)
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

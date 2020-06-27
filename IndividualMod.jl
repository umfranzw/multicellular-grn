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

    #genes, initial_protein_props = make_initial_genes(config)
    #genes = map(i -> GeneMod.rand_init(config, i, [ProteinPropsMod.Internal], [GeneMod.Id]), 1:config.run.num_initial_genes)

    genes = make_initial_genes(config)
    root_cell = Cell(config, genes)
    cell_tree = CellTree(root_cell)

    #initial_proteins = make_initial_proteins(config, initial_protein_props, length(genes), root_cell)
    initial_proteins = make_initial_proteins(config, genes, root_cell)

    indiv = Individual(config, genes, cell_tree, initial_proteins, 1.0, zeros(Int64, length(genes)))
    CellMod.insert_initial_proteins(root_cell, indiv.initial_cell_proteins)

    indiv
end

function make_initial_genes_preset(config::Config)
    genes = Array{Gene, 1}()
    initial_protein_props = Array{ProteinProps, 1}()
    genome_index = -1

    #ab
    ab_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ProteinPropsMod.Internal],
        tag=[UInt8(0)],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    ab_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Internal],
        tag=[UInt8(1)],
        fcn=ProteinPropsMod.Activate,
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )
    ab = pop_remaining_sites(config, genome_index, ab_bind_site, ab_prod_site)

    props = ProteinPropsMod.rand_init(
        config,
        type=[ab_bind_site.type],
        tag=[UInt8(0)],
        fcn=ProteinPropsMod.Activate,
    )
    push!(initial_protein_props, props)

    #bc
    bc_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ab_prod_site.type],
        tag=[ab_prod_site.tag],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    bc_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Application],
        tag=[UInt8(2)],
        action=[ProteinPropsMod.Divide],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )
    bc = pop_remaining_sites(config, genome_index, bc_bind_site, bc_prod_site)

    #ba
    ba_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ab_prod_site.type],
        tag=[ab_prod_site.tag],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    ba_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ab_bind_site.type],
        tag=[ab_bind_site.tag],
        fcn=ProteinPropsMod.Activate,
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )
    ba = pop_remaining_sites(config, genome_index, ba_bind_site, ba_prod_site)

    #de
    de_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ProteinPropsMod.Internal],
        tag=[UInt8(3)],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    plus_index = findfirst(i -> SymProbsMod.index_to_sym[i].val == :+, 1:length(SymProbsMod.index_to_sym)) - 1
    de_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Application],
        tag=[UInt8(4)],
        action=[ProteinPropsMod.SymProb],
        fcn=ProteinPropsMod.Activate,
        arg=[Int8(plus_index)],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )
    de = pop_remaining_sites(config, genome_index, de_bind_site, de_prod_site)

    props = ProteinPropsMod.rand_init(
        config,
        type=[de_bind_site.type],
        tag=[UInt8(3)],
        fcn=ProteinPropsMod.Activate
    )
    push!(initial_protein_props, props)

    #bf
    bf_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ab_prod_site.type],
        tag=[ab_prod_site.tag],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    bf_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Neighbour],
        tag=[UInt8(4)],
        fcn=ProteinPropsMod.Activate,
        arg=[Int8(4)], #right child
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )
    bf = pop_remaining_sites(config, genome_index, bf_bind_site, bf_prod_site)

    #bg
    bg_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ab_prod_site.type],
        tag=[ab_prod_site.tag],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    bg_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Neighbour],
        tag=[UInt8(5)],
        fcn=ProteinPropsMod.Activate,
        arg=[Int8(3)], #left child
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )
    bg = pop_remaining_sites(config, genome_index, bg_bind_site, bg_prod_site)

    #gh
    gh_bind_site = GeneMod.rand_bind_site(
        config,
        type=[bg_prod_site.type],
        tag=[bg_prod_site.tag],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )

    x_index = findfirst(i -> SymProbsMod.index_to_sym[i].val == :x, 1:length(SymProbsMod.index_to_sym)) - 1
    gh_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Application],
        tag=[UInt8(6)],
        action=[ProteinPropsMod.SymProb],
        fcn=ProteinPropsMod.Activate,
        arg=[Int8(x_index)],
        threshold=[IndividualMod.initial_threshold],
        consum_rate=[IndividualMod.initial_consum_rate]
    )
gh = pop_remaining_sites(config, genome_index, gh_bind_site, gh_prod_site)

#fnc
fnc_bind_site = GeneMod.rand_bind_site(
    config,
    type=[bf_prod_site.type], #neighbour-only receptor
    tag=[bf_prod_site.tag],
    threshold=[IndividualMod.initial_threshold],
    consum_rate=[IndividualMod.initial_consum_rate]
)

fnc_prod_site = GeneMod.rand_prod_site(
    config,
    type=[ProteinPropsMod.Neighbour],
    tag=[bc_prod_site.tag],
    fcn=ProteinPropsMod.Inhibit,
    arg=[Int8(1)], #parent
    threshold=[IndividualMod.initial_threshold],
    consum_rate=[IndividualMod.initial_consum_rate]
)
fnc = pop_remaining_sites(config, genome_index, fnc_bind_site, fnc_prod_site)

#fi
fi_bind_site = GeneMod.rand_bind_site(
    config,
    type=[bf_prod_site.type], #neighbour-only receptor
    tag=[bf_prod_site.tag],
    threshold=[IndividualMod.initial_threshold],
    consum_rate=[IndividualMod.initial_consum_rate]
)

one_index = findfirst(i -> SymProbsMod.index_to_sym[i].val == :1, 1:length(SymProbsMod.index_to_sym)) - 1
fi_prod_site = GeneMod.rand_prod_site(
    config,
    type=[ProteinPropsMod.Application],
    tag=[UInt8(7)],
    action=[ProteinPropsMod.SymProb],
    fcn=ProteinPropsMod.Activate,
    arg=[Int8(one_index)],
    threshold=[IndividualMod.initial_threshold],
    consum_rate=[IndividualMod.initial_consum_rate]
)
fi = pop_remaining_sites(config, genome_index, fi_bind_site, fi_prod_site)
push!(genes, gh, bg, ba, ab, bc, fnc, bf, fi, de)
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

function make_initial_genes(config::Config)
    genes = Array{Gene, 1}()

    for i in 1:config.run.num_initial_genes
        gene = GeneMod.rand_init(
            config,
            i,
            bind_site_types=[ProteinPropsMod.Internal],
            bind_logic=[GeneMod.Id]
        )
        push!(genes, gene)
    end

    genes
end

function make_initial_proteins(config::Config, genes::Array{Gene, 1}, root_cell::Cell)
    proteins = Array{Protein, 1}()
    num_concs = length(genes)
    half = (1.0 - IndividualMod.initial_threshold) / 2

    #choose genes that the initial proteins will be designed to bind to
    #note: these will only bind to bind sites (no inhibitory proteins here)
    indices = Random.shuffle(config.rng, 1:num_concs)
    for i in 1:config.run.num_initial_proteins
        index = indices[i]
        gene = genes[index]
        bind_site_index = RandUtilsMod.rand_int(config, 1, length(gene.bind_sites))
        bind_site = gene.bind_sites[bind_site_index]
        props = ProteinPropsMod.rand_init(
            config,
            type=[bind_site.type],
            tag=[bind_site.tag],
            fcn=ProteinPropsMod.Activate
        )
        protein = Protein(config, props, false, true, num_concs, root_cell.id)
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

function get_neighbour_at(cell::Cell, info::TreeInfo, loc::Int64)
    neighbour = nothing

    if loc == Int64(ProteinPropsMod.Top)
        if cell.parent != nothing
            neighbour = cell.parent
        end

    elseif loc == Int64(ProteinPropsMod.Left) || loc == Int64(ProteinPropsMod.Right)
        level = info.cell_to_level[cell]
        row = info.level_to_cell[level]
        cell_index = findall(c -> c == cell, row)[1]

        if loc == Int64(ProteinPropsMod.Left)
            if cell_index > 1
                neighbour = row[cell_index - 1]
            end
        else
            if cell_index < length(row)
                neighbour = row[cell_index + 1]
            end
        end

        #bottom (child)
    else
        num_children = length(cell.children)
        offset = loc - 3 + 1 #- 3 for left, top, right
        if num_children > 0 && offset <= num_children
            neighbour = cell.children[offset]
        end
    end

    neighbour
end

#returns cell's loc with respect to neighbour
function get_loc_relationship(cell::Cell, neighbour::Cell, neighbour_loc::Int64)
    rel = nothing

    if neighbour_loc == Int64(ProteinPropsMod.Top)
        offset = 1
        while offset <= length(neighbour.children) && cell != neighbour.children[offset]
            offset += 1
        end
        if offset <= length(neighbour.children)
            rel = Int8(2 + offset) #left=0, top=1, right=2, first_child=3, ...
        end

    elseif neighbour_loc == Int64(ProteinPropsMod.Left)
        rel = Int8(ProteinPropsMod.Right)

    elseif neighbour_loc == Int64(ProteinPropsMod.Right)
        rel = Int8(ProteinPropsMod.Left)

    else #bottom
        rel = Int8(ProteinPropsMod.Top)
    end

    rel
end

#transfers proteins to cell from neighbours
function run_neighbour_comm_for_cell(cell::Cell, info::TreeInfo)
    for src_loc in 0 : 2 + cell.config.run.max_children #the loc we're taking the proteins from
        neighbour = get_neighbour_at(cell, info, src_loc) #the cell we're taking the proteins from

        if neighbour != nothing
            #get neighbour's neighbour proteins with opposite locs (the ones that are being sent to this cell)
            #note: we use the id to ensure that we only get the ones that originated in the neighbour cell (it's possible they may have been transferred in from another neighbour)
            protein_loc = get_loc_relationship(cell, neighbour, src_loc)
            neighbour_proteins = ProteinStoreMod.get_neighbour_proteins_by_loc(neighbour.proteins, protein_loc, neighbour.id)

            if length(neighbour_proteins) > 0
                #compute the max amount of each protein that we can accept
                #accept_amount = cell.sensors[src_loc] / length(neighbour_proteins)
                for neighbour_protein in neighbour_proteins
                    #transfer_amount = min.(neighbour_protein.concs, accept_amount)
                    transfer_amount = neighbour_protein.concs

                    #add the neighbour protein into the source cell
                    dest_protein = ProteinStoreMod.get(cell.proteins, neighbour_protein.props)
                    if dest_protein == nothing
                        if ProteinStoreMod.num_proteins(cell.proteins) < cell.config.run.max_proteins_per_cell
                            dest_protein = Protein(cell.config, deepcopy(neighbour_protein.props), false, false, length(cell.gene_states), neighbour.id)
                            ProteinStoreMod.insert(cell.proteins, dest_protein)
                        end
                    end
                    if dest_protein != nothing
                        #cell.sensors[src_loc] -= transfer_amount
                        dest_protein.concs = clamp.(dest_protein.concs + transfer_amount, 0.0, 1.0)
                        neighbour_protein.concs -= transfer_amount
                    end

                    #note: if the neighbour_protein was completely tranferred, the decay step will delete it later
                end
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
            type = ProteinPropsMod.get_fcn(protein.props) == ProteinPropsMod.Inhibit ? ProdSite : BindSite
            for gs in cell.gene_states
                for i in 1:length(gs.bindings[type])
                    if gs.bindings[type][i] == protein
                        gs.bindings[type][i] = nothing
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
        prod_site.tag,
        prod_site.action,
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
            site = gene.bind_sites[site_index]
            eligible_proteins = get_bind_eligible_proteins_for_site(cell, gene, site)
            run_bind_for_site(indiv.config, cell.gene_states[gene_index], gene_index, BindSite, site_index, eligible_proteins)
        end

        for site_index in 1:length(gene.prod_sites)
            site = gene.prod_sites[site_index]
            eligible_proteins = get_bind_eligible_proteins_for_site(cell, gene, site)
            run_bind_for_site(indiv.config, cell.gene_states[gene_index], gene_index, ProdSite, site_index, eligible_proteins)
        end
    end
end

function get_bind_eligible_proteins_for_site(cell::Cell, gene::Gene, site::Union{BindSite, ProdSite})
    eligible_proteins = Array{Protein, 1}()
    search_proteins = ProteinStoreMod.get_by_types(cell.proteins, Set{ProteinPropsMod.ProteinType}( (ProteinPropsMod.Internal, ProteinPropsMod.Neighbour, ProteinPropsMod.Diffusion) ))
    
    if site isa BindSite
        #Binding logic for bind sites:
        # For proteins of type Internal:
        # -protein's type must match site's type
        # -protein's tag must match site's tag
        # -protein fcn must not be inhibit (inhibitory proteins bind to prod sites)
        # -protein's conc must be >= site's threshold

        #For proteins of type Neighbour:
        # -protein's type must match site's type
        # -protein's tag must match site's tag
        # -protein fcn must not be inhibit (inhibitory proteins bind to prod sites)
        # -protein's conc must be >= site's threshold
        # -protein's src_cell_id must not point to the current cell (to prevent self-binding)

        #For proteins of type Diffusion:
        # -protein's type must match site's type
        # -protein's tag must match site's tag
        # -protein fcn must not be inhibit (inhibitory proteins bind to prod sites)
        # -protein's conc must be >= site's threshold
        # -protein's src_cell_id must not point to the current cell (to prevent self-binding)

        #Proteins of type Application do not bind
        
        for protein in search_proteins
            #check the conditions that all protein types have in common
            #note: this completes all the checks needed for proteins of type internal
            eligible = (protein.props.type == site.type &&
                        protein.props.tag == site.tag &&
                        ProteinPropsMod.get_fcn(protein.props) != ProteinPropsMod.Inhibit &&
                        protein.concs[gene.genome_index] >= site.threshold)

            #check conditions specific to neighbour and diffusion types
            if eligible
                #note: leave these two cases separate for now in case they diverge later...
                if protein.props.type == ProteinPropsMod.Neighbour
                    eligible = protein.src_cell_id != cell.id
                    
                elseif protein.props.type == ProteinPropsMod.Diffusion
                    eligible = protein.src_cell_id != cell.id
                end
            end

            if eligible
                push!(eligible_proteins, protein)
            end
        end

    else
        #Binding logic for prod sites (inhibitory):
        # For proteins of type Internal:
        # -protein's tag must match site's tag
        # -protein fcn must be inhibit (only inhibitory proteins bind to prod sites)
        # -site's fcn must not be inhibitory (inhibitory proteins can't block the production of inhibitory proteins)
        # -protein's conc must be >= site's threshold

        #For proteins of type Neighbour:
        # -protein's tag must match site's tag
        # -protein fcn must be inhibit (only inhibitory proteins bind to prod sites)
        # -site's fcn must not be inhibitory (inhibitory proteins can't block the production of inhibitory proteins)
        #   -this also prevents neighbour proteins from binding to same column they were produced in
        # -protein's src_cell_id must not point to the current cell (to prevent self-binding)

        #For proteins of type Diffusion:
        # -protein's tag must match site's tag
        # -protein fcn must be inhibit (only inhibitory proteins bind to prod sites)
        # -site's fcn must not be inhibitory (inhibitory proteins can't block the production of inhibitory proteins)
        #   -this also prevents diffusion proteins from binding to same column they were produced in
        # -protein's src_cell_id must not point to the current cell (to prevent self-binding)

        #Proteins of type Application do not bind
        
        for protein in search_proteins
            #check the conditions that all protein types have in common
            #note: this completes all the checks needed for proteins of type internal
            eligible = (protein.props.tag == site.tag &&
                        ProteinPropsMod.get_fcn(protein.props) == ProteinPropsMod.Inhibit &&
                        ProteinPropsMod.get_fcn(site.arg) != ProteinPropsMod.Inhibit &&
                        protein.concs[gene.genome_index] >= site.threshold)
            
            #check conditions specific to neighbour and diffusion types
            if eligible
                #note: leave these two cases separate for now in case they diverge later...
                if protein.props.type == ProteinPropsMod.Neighbour
                    eligible = protein.src_cell_id != cell.id
                    
                elseif protein.props.type == ProteinPropsMod.Diffusion
                    eligible = protein.src_cell_id != cell.id
                end
            end

            if eligible
                push!(eligible_proteins, protein)
            end
        end
    end
    
    eligible_proteins
end

function run_bind_for_site(config::Config, gs::GeneState, gene_index::Int64, site_type::Union{Type{BindSite}, Type{ProdSite}}, site_index::Int64, eligible_proteins::Array{Protein, 1})
    run_bind_for_site_max(config, gs, gene_index, site_type, site_index, eligible_proteins)
end

function run_bind_for_site_max(config::Config, gs::GeneState, gene_index::Int64, site_type::Union{Type{BindSite}, Type{ProdSite}}, site_index::Int64, eligible_proteins::Array{Protein, 1})
    if length(eligible_proteins) > 0
        max_protein = nothing
        for protein in eligible_proteins
            protein_conc = protein.concs[gene_index]
            if max_protein == nothing || protein_conc > max_protein.concs[gene_index]
                max_protein = protein
            end
        end

        GeneStateMod.bind(gs, max_protein, site_type, site_index)
        
    else
        state = GeneStateMod.get_binding(gs, site_type, site_index)
        if state != nothing
            GeneStateMod.unbind(gs, site_type, site_index)
        end
    end
end

function run_bind_for_site_roulette(config::Config, gs::GeneState, gene_index::Int64, site_type::Union{Type{BindSite}, Type{ProdSite}}, site_index::Int64, eligible_proteins::Array{Protein, 1})
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
        GeneStateMod.bind(gs, sel_protein, site_type, site_index)

    else
        state = GeneStateMod.get_binding(gs, site_type, site_index)
        if state != nothing
            #@info @sprintf("Protein unbinding")
            GeneStateMod.unbind(gs, site_type, site_index)
        end
    end
end

end

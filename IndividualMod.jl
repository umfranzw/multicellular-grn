module IndividualMod

using RunMod
using GeneMod
using GeneStateMod
using ProteinMod
using ProteinPropsMod
using CellMod
using ProteinStoreMod
using SymMod
using DiffusionMod
using CellTreeMod

import Random
import RandUtilsMod
import Base.show

export Individual,
    rand_init, run_bind

struct Individual
    run::Run
    genes::Array{Gene, 1}
    root_cell::Cell
    initial_cell_proteins::Array{Protein, 1}
    #note: this is a value in [0.0, 1.0], where 0.0 is optimal
    fitness::Float64
end

function show(io::IO, indiv::Individual)
    println(io, "This is an individual.")
end

function rand_init(run::Run)
    genes = map(i -> GeneMod.rand_init(run, i), 1:run.num_genes)
    root_cell = Cell(run, genes, nothing, Sym(:x, SymMod.DataVar))
    
    initial_proteins = Array{Protein, 1}()
    
    for i in 1:run.num_initial_proteins
        #all initial proteins will be type Reg, and have target Intra
        #this forces inter-cell communication and application proteins to be produced by a network
        type = ProteinPropsMod.Reg
        target = ProteinPropsMod.Intra
        reg_action = RandUtilsMod.rand_enum_val(run, ProteinPropsMod.ProteinRegAction)
        app_action = RandUtilsMod.rand_enum_val(run, ProteinPropsMod.ProteinAppAction)
        
        protein = Protein(run, ProteinProps(type, target, reg_action, app_action), true)
        push!(initial_proteins, protein)

        #it is possible that not all initial proteins in the array are unique. That's ok, since they'll be subject to evolution.
        #However, we need to ensure that we only insert unique proteins into the root cell's store.
        #note: this logic means some individual's root cells may have fewer initial proteins than others...
        if !ProteinStoreMod.contains(root_cell.proteins, protein)
            #note: we push a copy so the indiv's initial_cell_proteins array stays intact as the simulation modifies protein's concs
            #in the root cell
            ProteinStoreMod.insert(root_cell.proteins, ProteinMod.copy(protein), false)
        end
    end
    
    Individual(run, genes, root_cell, initial_proteins, 1.0)
end

#resets everything to the way it was before the reg sim (so
#we can run the reg sim again on the next ea_step)
function reset(indiv::Individual)
    #just re-initialize the cell (this discards the rest of the tree, along with any protein bindings)
    indiv.root_cell = Cell(run, genes, nothing, Sym(:x, SymMod.DataVar))
    #re-insert (copies of) the initial proteins
    for protein in indiv.initial_proteins
        if !ProteinStoreMod.contains(root_cell.proteins, protein)
            ProteinStoreMod.insert(root_cell.proteins, ProteinMod.copy(protein), false)
        end
    end
end

function run_protein_app(indiv::Individual)
    #we want to visit cells in breadth-first order
    #build an array of the cells (in bfs order) so that as the tree is modified,
    #we don't get messed up by any modifications (eg. new nodes that get added)
    bfs_list = Array{Cell}()
    CellTree.bf_traverse(indiv.root_cell, c -> push!(bfs_list, c))

    for cell in bfs_list
        run_app_for_cell(cell)
    end
end

function run_protein_app_for_cell(cell::Cell)
    #get all proteins (from this cell) that are eligible for application
    app_proteions = ProteinStoreMod.get_by_type(cell.proteins, ProteinPropsMod.App)

    #build a list of tuples of the form (protein, sum of concs), where each protein has a sum >= protein_app_threshold
    pairs = Array{Tuple{Protein, Float64}}()
    for protein in app_proteins
        conc_sum = sum(protein.concs)
        if conc_sum >= cell.run.protein_app_threshold
            push!(pairs, (p, conc_sum))
        end
    end
    #sort in descending order by sum - we'll apply they in this order
    sort!(pairs; by=p -> p[2], rev=true)

    #apply the proteins
    for pair in pairs
        protein = pair[1]
        ProteinAppActionMod.run_app_action(cell, protein)
    end
end

function run_bind(indiv::Individual)
    CellTreeMod.traverse(indiv.root_cell, cell -> run_bind_for_cell(indiv, cell))
end

function run_produce(indiv::Individual)
    CellTreeMod.traverse(indiv.root_cell, cell -> run_produce_for_cell(indiv, cell))
end

function run_diffuse(indiv::Individual)
    DiffusionMod.diffuse_intra_cell_proteins(indiv.root_cell)
    DiffusionMod.diffuse_inter_cell_proteins(indiv.root_cell)
end

function run_decay(indiv::Individual)
    CellTreeMod.traverse(indiv.root_cell, cell -> run_decay_for_cell(cell))
end

function is_decayed(protein::Protein, thresh::Float64)
    i = 1
    while protein.concs[i] < thresh && i <= length(protein.concs)
        i += 1
    end

    #if the loop index made it past the length of the protein, all the concs are below the threshold (i.e. the protein is decayed)
    i > length(protein.concs)
end

function run_decay_for_cell(cell::Cell)
    proteins = ProteinStoreMod.get_all(cell.proteins)
    for protein in proteins
        #decrease the concentration using decay_rate
        protein.concs = max.(protein.concs - protein.concs * cell.run.decay_rate, zeros(length(protein.concs)))

        #remove any proteins that have decayed below the allowable threshold
        if is_decayed(protein, cell.run.min_protein_conc)
            #note: as long as the decay_threshold < reg_decay_threshold, and run_bind() executes before run_decay(),
            #no decayed protein will be bound to anything at this point, so we don't need to go through the gene
            #states to remove bindings
            cell.proteins.remove(protein)
        end
    end
end

function run_produce_for_cell(indiv::Individual, cell::Cell)
    for gene_index in 1:indiv.run.num_genes
        gene_state = cell.gene_states[gene_index]
        gene = indiv.genes[gene_index]
        rates = GeneStateMod.get_prod_weights(gs)

        #For intra prod site
        run_produce_for_site(cell, gene_index, GeneMod.Intra, rates.intra)
        #For inter prod site
        run_produce_for_site(cell, gene_index, GeneMod.Inter, rates.inter)
    end
end

function run_produce_for_site(cell::Cell, gene_index::Int64, site_type::GeneMod.ProdSites, rate::Float64)
    #get the props for the protein that will be produced
    props = gene.prod_sites[site_type]
    #check if protein already exists in this cell's store
    protein = ProteinStoreMod.get(cell.proteins, props)
    #if not, create and insert it
    if protein == nothing
        #note: protein will be initialized with conc values of zero
        protein = Protein(run, props, false)
        ProteinStoreMod.insert(cell.proteins, protein, true)
    end
    
    #increment the conc using the rate
    #note: this will only increment the conc directly over the gene
    #the diffusion will spread this out later
    protein.concs[gene_index] = min(protein.concs[gene_index] + rate, 1.0)
end

function run_bind_for_cell(indiv::Individual, cell::Cell)
    for gene_index in 1:indiv.run.num_genes
        gene = indiv.genes[gene_index]
        
        #run binding for each of the regulatory sites
        for site_index in 1:length(gene.reg_sites)
            site = gene.reg_sites[site_index]
            site_type = GeneMod.RegSites(site_index)

            #Get eligable proteins
            #for site types GeneMod.IntraIntra and GeneMod.IntraInter
            if site.target == ProteinPropsMod.Intra
                eligable_proteins = get_bind_eligable_proteins_for_intra_site(cell.proteins, gene_index, site, indiv.run.reg_bind_threshold)
                
            #for site types GeneMod.InterIntra and GeneMod.InterInter
            else
                proteins = get_bind_eligable_proteins_for_inter_site(cell.proteins, gene_index, site, indiv.run.reg_bind_threshold)
            end
            
            #run the binding
            run_bind_for_site(cell.gene_states[gene_index], site_type, eligable_proteins)
        end
    end
end

function get_bind_eligable_proteins_for_intra_site(ps::ProteinStore, gene_index::Int64, site::ProteinProps, bind_threshold::Float64)
    proteins = values(ProteinStore.get_by_target(ps, ProteinPropsMod.Intra))

    #for props:
    #- need to ensure that protein's type matches site type (Reg)
    #- protein's target will match site target (Intra) due the the filtering from the call above
    #- protein's reg action doesn't need to match (could be Activate or Inhibit)
    #- protein's app action is irrelevant (since this is a reg protein & reg site)
    filter(p -> p.concs[gene_index] >= bind_threshold && p.props.type == site.type, proteins)
end

function get_bind_eligable_proteins_for_inter_site(ps::ProteinStore, gene_index::Int64, site::ProteinProps, bind_threshold::Float64)
    inter_cell_proteins = values(ProteinStore.get_by_target(ps, ProteinPropsMod.Inter))
    eligable_proteins = Array{Protein, 1}()

    for protein in inter_cell_proteins
        #we want to make sure that inter-cell proteins don't bind to genes in the cell that produced them (don't "self-bind")
        #This call checks if the given store owns the given protein.
        #note: it's possible that the store owns the protein, AND the protein was ALSO produced by another cell.
        #Since we can't differentiate between the two cases, we still don't allow the protein to bind if both are true.
        is_self_binding = ProteinStore.is_owned_intercell_protein(ps, protein)
        above_thresh = protein.concs[gene_index] >= bind_threshold

        #for props:
        #- need to ensure that protein's type matches site type (Reg)
        #- protein's target will match site target (Inter) due the the filtering from the call above
        #- protein's reg action doesn't need to match (could be Activate or Inhibit)
        #- protein's app action is irrelevant (since this is a reg protein & reg site)
        seq_matches = protein.props.type == site.type

        if !is_self_binding && above_thresh && seq_matches
            push!(eligable_proteins, protein)
        end
    end

    eligable_proteins
end

function run_bind_for_site(gs::GeneState, site_type::Union{GeneMod.RegSites, GeneMod.ProdSites}, eligable_proteins::Array{Protein, 1})
    if length(eligible_proteins) > 0
        #use roulette wheel style selection to pick the protein
        conc_sum = foldl((s, p) -> s + p.concs[gene_index], eligible_proteins; init=0.0)
        next_bound = 0.0
        wheel = []
        for protein in eligible_proteins
            next_bound = protein.concs[gene_index] / conc_sum
            push!(wheel, next_bound)
        end

        wheel[end] = 1.0 #just in case we've got some floating point error

        sel_index = 1
        r = RandUtilsMod.rand_float(run) #random value in [0, 1)
        while r >= wheel[sel_index]
            sel_index += 1
        end

        sel_protein = eligible_proteins[sel_index]

        GeneState.bind(gene_state, sel_protein, site_type)

    else
        GeneState.unbind(gene_state, site_type)
    end
end

end

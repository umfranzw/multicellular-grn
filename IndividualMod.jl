module IndividualMod

using RunMod
using GeneMod
using GeneStateMod
using ProteinMod
using CellMod
using ProteinStoreMod
using SymMod
using DiffusionMod

import Random
import RandUtilsMod

export Individual,
    rand_init, run_bind

struct Individual
    run::Run
    genes::Array{Gene, 1}
    cells::Array{Cell, 1}
    initial_cell_proteins::Array{Protein, 1}
end

function rand_init(run::Run)
    genes = map(i -> GeneMod.rand_init(run, i), 1:run.num_genes)
    initial_cell = Cell(run, genes, nothing, Sym(:x, SymMod.DataVar))
    
    initial_proteins = Array{Protein, 1}()
    
    for i in 1:run.num_initial_proteins
        type = RandUtilsMod.rand_enum_val(ProteinMod.ProteinType)
        target = RandUtilsMod.rand_enum_val(ProteinMod.ProteinTarget)
        reg_action = RandUtilsMod.rand_enum_val(ProteinMod.ProteinRegAction)
        growth_action = RandUtilsMod.rand_enum_val(ProteinMod.ProteinGrowthAction)
        
        protein = Protein(run, ProteinProps(type, target, reg_action, growth_action), true)
        
        #make sure the initial proteins are all unique
        #note: this logic means some cells in the population may have fewer initial proteins than others...
        if !ProteinStoreMod.contains(protein)
            push!(initial_proteins, protein)
            #note: we push a copy so the indiv's array stays intact as the simulation modifies protein's concs
            ProteinStoreMod.insert(initial_cell.proteins, ProteinMod.copy(protein), false)
        end
    end
    
    Individual(run, genes, [initial_cell], initial_proteins, store)
end

# function divide_cell(indiv::Individual, src_cell::Cell, dest_index::Int64)
#     dest_cell = Cell(indiv.run, indiv.genes)

#     ProteinStoreMod.insert_cell(indiv.protein_store, dest_index)
#     insert!(indiv.cell, dest_index)

#     #insert the initial proteins for the new cell
#     #src cell donates half its intra-cell proteins' concentrations to the dest cell
#     intra_proteins = ProteinStoreMod.get_proteins_by_scope(indiv.protein_store, ProteinMod.IntraCell)

#     for protein in intra_proteins
#         protein.concs[src_cell] /= 2
#         protein.concs[dest_cell] = copy(protein.concs[src_cell])
#     end

#     src_cell.energy /= 2
# end

# function run_growth(indiv::Individual)
#     #as we insert, the list will grow - we need to compensate for this
#     cells = copy(indiv.cells) #we will loop through this copy so that we don't run into cells inserted on a previous iter
#     index_offset = 0 #offset for cells already inserted

#     for i in 1:length(cells)
#         prob, dir = calc_cell_divide_prob(cells[i])
#         if RandUtils.rand_float(indiv.run) < prob
#             if dir == ProteinMod.GrowUp
#                 divide_cell(indiv, cells[i], max(i + index_offset - 1, 1))
#             elseif dir == ProteinMod.GrowDown
#                 divide_cell(indiv, cells[i], min(i + index_offset + 1, length(indiv.cells)))
#             end
#             index_offset += 1
#         end
#     end
# end

function run_bind(indiv::Individual)
    for i in 1:length(indiv.cells)
        run_bind_for_cell(indiv, i)
    end
end

function run_produce(indiv::Individual)
    for i in 1:length(indiv.cells)
        run_produce_for_cell(indiv, i)
    end
end

function run_diffuse(indiv::Individual)
    DiffusionMod.diffuse_intra_cell_proteins(indiv)
    DiffusionMod.diffuse_inter_cell_proteins(indiv)
end

function run_produce_for_cell(indiv::Individual, cell_index::Int64)
    cell = indiv.cells[cell_index]
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
    protein = ProteinStoreMod.get(cell.protein_store, props)
    #if not, create and insert it
    if protein == nothing
        #note: protein will be initialized with conc values of zero
        protein = Protein(run, props, false)
        ProteinStoreMod.insert(cell.protein_store, protein, true)
    end
    
    #increment the conc using the rate
    #note: this will only increment the conc directly over the gene
    #the diffusion will spread this out later
    protein.concs[gene_index] = min(protein.concs[gene_index] + rate, 1.0)
end

function run_bind_for_cell(indiv::Individual, cell_index::Int64)
    cell = indiv.cells[cell_index]
    for gene_index in 1:indiv.run.num_genes
        gene = indiv.genes[gene_index]
        
        #run binding for each of the regulatory sites
        for site_index in 1:length(gene.reg_sites)
            site = gene.reg_sites[site_index]
            site_type = GeneMod.RegSites(site_index)

            #Get eligable proteins
            #for site types GeneMod.IntraIntra and GeneMod.IntraInter
            if site.target == ProteinMod.Intra
                eligable_proteins = get_bind_eligable_proteins_for_intra_site(cell.protein_store, gene_index, site, indiv.run.reg_bind_threshold)
                
            #for site types GeneMod.InterIntra and GeneMod.InterInter
            else
                proteins = get_bind_eligable_proteins_for_inter_site(cell.protein_store, cell_index, gene_index, site, indiv.run.reg_bind_threshold)
            end
            
            #run the binding
            run_bind_for_site(cell.gene_states[gene_index], site_type, eligable_proteins)
        end
    end
end

function get_bind_eligable_proteins_for_intra_site(ps::ProteinStore, gene_index::Int64, site::ProteinProps, bind_threshold::Float64)
    proteins = values(ProteinStore.get_by_target(ps, ProteinMod.Intra))

    #for props:
    #- need to ensure that protein's type matches site type (Reg)
    #- protein's target will match site target (Intra) due the the filtering from the call above
    #- protein's reg action doesn't need to match (could be Activate or Inhibit)
    #- protein's growth action is irrelevant (since this is a reg protein & reg site)
    filter(p -> p.concs[gene_index] >= bind_threshold && p.props.type == site.type, proteins)
end

function get_bind_eligable_proteins_for_inter_site(ps::ProteinStore, gene_index::Int64, site::ProteinProps, bind_threshold::Float64)
    inter_cell_proteins = values(ProteinStore.get_by_target(ps, ProteinMod.Inter))
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
        #- protein's growth action is irrelevant (since this is a reg protein & reg site)
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

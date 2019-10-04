module IndividualMod

using RunMod
using GeneMod
using GeneStateMod
using ProteinMod
using CellMod
using ProteinStoreMod
using SymMod

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

function run_regulate(indiv::Individual)
    for i in 1:length(indiv.cells)
        run_regulate_for_cell(indiv, i)
    end
end

function run_diffuse(indiv::Individual)
    DiffusionMod.diffuse_intra_cell_proteins(indiv)
    DiffusionMod.diffuse_inter_cell_proteins(indiv)
end

function run_regulate_for_cell(indiv::Individual, cell_index::Int64)
    cell = indiv.cells[cell_index]
    for gene_index in 1:indiv.run.num_genes
        gene_state = cell.gene_states[gene_index]

        bound_protein = gene_state.reg_site_binding
        if bound_protein != nothing
            affinity = ProteinMod.get_bind_affinity(bound_protein)
            if affinity == ProteinMod.Reg
                reg_action = ProteinMod.get_reg_action(bound_protein)
                if reg_action == ProteinMod.RateUp
                    gene_state.prod_rate = min(gene_state.prod_rate + indiv.run.prod_rate_incr, 1.0)
                elseif reg_action == ProteinMod.RateDown
                    gene_state.prod_rate = max(gene_state.prod_rate - indiv.run.prod_rate_incr, 0.0)
                end
            end
        end
    end
end

function run_produce_for_cell(indiv::Individual, cell_index::Int64)
    cell = indiv.cells[cell_index]
    for gene_index in 1:indiv.run.num_genes
        gene_state = cell.gene_states[gene_index]
        gene = indiv.genes[gene_index]
        
        for site_index in 1:indiv.run.num_bind_sites
            #check if something's bound to the bind site and nothing's bound to the prod site
            bound_protein = gene_state.bind_site_bindings[site_index]
            inhibit_protein = gene_state.prod_site_bindings[site_index]
            if bound_protein != nothing && inhibit_protein == nothing
                product_seq = gene.bind_sites[site_index]
                
                #if product is not yet in store, create it
                product_protein = ProteinStoreMod.get_protein(indiv.protein_store, product_seq)
                if product_protein == nothing
                    product_protein = ProteinMod.Protein(run, product_seq, length(indiv.cells), false)
                    ProteinStoreMod.insert_protein(indiv.protein_store, product_protein, cell)

                #otherwise, make sure this cell is a src
                else
                    #note: don't need to worry about inserting duplicates here, since these are stored in a set
                    ProteinStoreMod.add_src_cell(indiv.protein_store, cell)
                end

                #increment conc above this gene (diffusion will spread it out later)
                product_protein.concs[cell_index][gene_index] = min(produce_protein.concs[cell_index][gene_index] + gene_state.output_rate, 1.0)
            end
        end
    end
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
            
            #do binding
            run_bind_for_site(cell.gene_states[gene_index], site_type, eligable_proteins)
        end
    end
end

function get_bind_eligable_proteins_for_intra_site(ps::ProteinStore, gene_index::Int64, site::ProteinProps, bind_threshold::Float64)
    proteins = values(ProteinStore.get_by_target(ps, ProteinMod.Intra))

    filter(p -> p.concs[gene_index] >= bind_threshold && p.props == site, proteins)
end

function get_bind_eligable_proteins_for_intra_site(ps::ProteinStore, gene_index::Int64, site::ProteinProps, bind_threshold::Float64)
    inter_cell_proteins = values(ProteinStore.get_by_target(ps, ProteinMod.Inter))
    eligable_proteins = Array{Protein, 1}()

    for protein in inter_cell_proteins
        #we need to make sure that inter-cell proteins don't bind to genes in the cell that produced them (don't "self-bind")
        is_self_binding = ProteinStore.is_owned_intercell_protein(ps, protein)
        above_thresh = protein.concs[gene_index] >= bind_threshold
        seq_matches = protein.props == site

        if !is_self_binding && above_thresh && seq_matches
            push!(eligable_proteins, protein)
        end
    end

    eligable_proteins
end

function run_bind_for_site(gs::GeneState, site_type::Union{GeneMod.RegSites, GeneMod.ProdSites}, eligable_proteins::Array{Protein, 1})
    if length(eligible_proteins) > 0
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

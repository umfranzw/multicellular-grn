module BindStepMod

using RunMod
using IndividualMod
using CellTreeMod
using CellMod
using GeneMod
using ProteinStoreMod
using ProteinPropsMod
using ProteinMod
using GeneStateMod
using RandUtilsMod

function run_bind(indiv::Individual)
    CellTreeMod.traverse(cell -> run_bind_for_cell(indiv, cell), indiv.cell_tree)
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

module GrowthMod

using RunMod
using IndividualMod
using GeneMod
using ProteinMod

import RandUtilsMod
import Random
import GenomeUtilsMod
import RegSimInfoMod
import ProteinPropsMod
import CellMod

@enum LinkMethod::Int8 FromInitialProtein FromProdSite

function grow(run::Run, pop::Array{Individual, 1})
    if run.multithreaded
        Threads.@threads for indiv in pop
            if can_grow(indiv)
                grow_indiv(indiv)
            end
        end
    else
        for indiv in pop
            if can_grow(indiv)
                grow_indiv(indiv)
            end
        end
    end
end

function can_grow(indiv::Individual)
    length(indiv.genes) < indiv.config.run.max_genes &&
        RegSimInfoMod.get_bind_coverage(indiv.reg_sim_info) >= indiv.config.run.growth_threshold &&
        RandUtilsMod.rand_decision(indiv.config, indiv.config.run.growth_prob)
end

function grow_indiv(indiv::Individual)
    num_concs = length(indiv.genes)
    half = (1.0 - indiv.config.run.bind_threshold) / 2
    
    #add a new initial protein, if possible
    if length(indiv.initial_cell_proteins) < indiv.config.run.max_initial_proteins
        props = ProteinPropsMod.rand_init(
            indiv.config,
            type=[ProteinPropsMod.Internal],
            fcn=ProteinPropsMod.Activate #TODO: allow initial proteins (aside from the very first one) to be inhibitory
        )
        protein = Protein(indiv.config, props, false, true, num_concs, indiv.cell_tree.root.id)
        protein.concs = RandUtilsMod.rand_floats(indiv.config, indiv.config.run.bind_threshold + half, 1.0, num_concs)
        push!(indiv.initial_cell_proteins, protein)
        CellMod.insert_initial_proteins(indiv.cell_tree.root, [protein])
    end

    #determine the index of the new gene
    gene_index = RandUtilsMod.rand_int(indiv.config, 1, length(indiv.genes) + 1)

    #determine how the new gene will be linked into the genome
    link_method = RandUtilsMod.rand_enum_val(indiv.config, LinkMethod)

    #get the properties we'll need to link against
    link_info = nothing
    if link_method == FromInitialProtein
        #choose an initial protein and get its props
        link_info = Random.rand(indiv.config.rng, indiv.initial_cell_proteins).props
    elseif link_method == FromProdSite
        #choose a prod site and get its props
        prod_sites = Array{ProdSite, 1}()
        for gene in indiv.genes
            for prod_site in gene.prod_sites
                #can't link against sites that produce application proteins
                if prod_site.type != ProteinPropsMod.Application
                    push!(prod_sites, prod_site)
                end
            end
        end

        if length(prod_sites) > 0
            link_info = Random.rand(indiv.config.rng, prod_sites)
        end
    end
    
    #build the gene
    if link_info != nothing
        gene = GeneMod.rand_init(
            indiv.config,
            gene_index,
            #bind site type must match the type of the things we're linking with
            bind_site_types=[link_info.type],
            bind_logic=[GeneMod.Id]
        )

        GenomeUtilsMod.insert_genes(indiv, [gene], gene_index)
    end
end


end

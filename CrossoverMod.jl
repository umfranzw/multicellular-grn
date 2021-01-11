module CrossoverMod

#Assumes selection has been done beforehand (pop contains the newly selected individuals)

using RunMod
using IndividualMod
using GeneMod
using ProteinMod
import ProteinPropsMod
import Random

config = nothing

function crossover(run::Run, pop::Array{Individual, 1}, ea_step::Int64)
    global config
    
    if config == nothing
        config = Config(run, UInt64(length(pop) + 1))
    end

    #use crossover_proportion param (eg. 40% of pop): pop_size * 0.4 rounded to the nearest even number
    num_cross = Int64(floor(run.pop_size * run.cross_prop))
    num_pairs = num_cross ÷ 2
    #note: duplication of fitter individuals has already occurred in the selection
    #So here, we only allow each individual to participate in crossover once
    #This also ensures that children are not crossed on iterations after the one they were generated on
    p1_choices = Random.shuffle(config.rng, 1:length(pop))
    p2_choices = Random.shuffle(config.rng, 1:length(pop))
    for i in 1:num_pairs
        p1_index = pop!(p1_choices)
        p2_index = pop!(p2_choices)
        p1 = pop[p1_index]
        p2 = pop[p2_index]

        pop_seed_offset = run.pop_size + 1
        ea_step_seed_offset = ea_step * num_cross
        p1_to_p2_ranges = pick_ranges(p1, p2)
        c1 = nothing

        if p1_to_p2_ranges != nothing
            p1_range, p2_range = p1_to_p2_ranges
            c1_offset = (i - 1) * 2 + 1
            c1 = make_child(p1, p2, p1_range, p2_range, pop_seed_offset + ea_step_seed_offset + c1_offset)
        end

        p2_to_p1_ranges = pick_ranges(p2, p1)
        c2 = nothing
        if p2_to_p1_ranges != nothing #if there's a way to link p2 to p1
            p2_range, p1_range = p2_to_p1_ranges
            c2_offset = i * 2 #note: this is equivalent to (i - 1) * 2 + 2
            c2 = make_child(p2, p1, p2_range, p1_range, pop_seed_offset + ea_step_seed_offset + c2_offset)
        end

        # if i == 1 && (c1 != nothing || c2 != nothing)
        #     println("p1:")
        #     println(IndividualMod.get_gene_descs_str(p1))
        #     println("p2:")
        #     println(IndividualMod.get_gene_descs_str(p2))
        # end
        
        if c1 != nothing
            pop[p1_index] = c1
            # if i == 1
            #     println("c1:")
            #     println(p1_to_p2_ranges)
            #     println(IndividualMod.get_gene_descs_str(c1))
            # end
        end
        if c2 != nothing
            pop[p2_index] = c2
            # if i == 1
            #     println("c2:")
            #     println(p2_to_p1_ranges)
            #     println(IndividualMod.get_gene_descs_str(c2))
            # end
        end
    end
end

function build_site_dicts(genes::Array{Gene, 1})
    #(type, tag) => array of genome_indices
    prod_sites = Dict{Tuple{ProteinPropsMod.ProteinType, UInt8}, Array{Int64, 1}}()
    bind_sites = Dict{Tuple{ProteinPropsMod.ProteinType, UInt8}, Array{Int64, 1}}()
    for gene in genes
        for prod_site in gene.prod_sites
            key = (prod_site.type, prod_site.tag)
            if key ∉ keys(prod_sites)
                prod_sites[key] = Array{Int64, 1}()
            end
            push!(prod_sites[key], gene.genome_index)
        end
        
        for bind_site in gene.bind_sites
            key = (bind_site.type, bind_site.tag)
            if key ∉ keys(bind_sites)
                bind_sites[key] = Array{Int64, 1}()
            end
            push!(bind_sites[key], gene.genome_index)
        end
    end

    bind_sites, prod_sites
end

function pick_ranges(p1::Individual, p2::Individual)
    ranges = nothing
    pts = find_link_points(p1, p2)
    if length(pts) > 0
        p1_end, p2_start = Random.rand(p1.config.rng, pts)
        ranges = (1:p1_end, p2_start:length(p2.genes))
    end

    ranges
end

#find all the ways we can link p1 to p2
function find_link_points(p1::Individual, p2::Individual)
    links = Array{Tuple{Int64, Int64}, 1}() #[(p1_genome_index, ps_genome_index), ...]

    #build dictionaries to make lookups fast & easy later on
    #(type, tag) => array of genome_indices
    p1_bind_sites, p1_prod_sites = build_site_dicts(p1.genes)
    p2_bind_sites, p2_prod_sites = build_site_dicts(p2.genes)

    function append_links(
        src_key::Tuple{ProteinPropsMod.ProteinType, UInt8}, #(type, tag)
        src_sites::Dict{Tuple{ProteinPropsMod.ProteinType, UInt8}, Array{Int64, 1}},
        dest_sites::Dict{Tuple{ProteinPropsMod.ProteinType, UInt8}, Array{Int64, 1}}
    )
        if src_key in keys(dest_sites) #src_key is (protein_type, tag)
            #can link every src gene that produces the src_key protein with every dest gene that accepts it
            for src_index in src_sites[src_key]
                for dest_index in dest_sites[src_key]
                    push!(links, (src_index, dest_index))
                end
            end
        end
    end
    
    links = Array{Tuple{Int64, Int64}, 1}(); #[(p1_gene_index, p2_gene_index), ...]
    for prod_key in keys(p1_prod_sites)
        #p1_prod_site -> p2_bind_site
        append_links(prod_key, p1_prod_sites, p2_bind_sites)

        #p1_prod_site -> p2_prod_site (for inhibitory proteins)
        append_links(prod_key, p1_prod_sites, p2_prod_sites)
    end

    links
end

function make_child(p1::Individual, p2::Individual, p1_gene_range::UnitRange{Int64}, p2_gene_range::UnitRange{Int64}, seed_offset::Int64)
    #note: although both children are replaced in the current population, it's possible that one of them may be rolled back to
    #its parent at the end of this ea_step (see EvAlgMod.enforce_fitness_front()). So we can't just steal the genes and config
    #from the parents - we have to deepcopy them.
    genes = Array{Gene, 1}()
    initial_proteins = Array{Protein, 1}()
    genome_index = 1
    total_concs = (p1_gene_range.stop - p1_gene_range.start + 1) + (p2_gene_range.stop - p2_gene_range.start + 1) # == total_genes
    for i in p1_gene_range
        gene = p1.genes[i]
        child_gene = deepcopy(gene)
        child_gene.genome_index = genome_index
        push!(genes, child_gene)
        genome_index += 1

        #add a matching initial protein, if possible
        if i <= length(p1.initial_cell_proteins) && length(initial_proteins) < p1.config.run.max_initial_proteins
            new_protein = deepcopy(p1.initial_cell_proteins[i])
            resize!(new_protein.concs, total_concs)
            #make sure no concs overlap with p2 genes
            new_protein.concs[p1_gene_range.stop + 1:total_concs] .= 0 #note: it's ok if the start point is outside the array - julia doesn't care
            push!(initial_proteins, new_protein)
        end
    end

    for i in p2_gene_range
        gene = p2.genes[i]
        child_gene = deepcopy(gene)
        child_gene.genome_index = genome_index
        push!(genes, child_gene)
        genome_index += 1

        #add a matching initial protein, if possible
        if i <= length(p2.initial_cell_proteins) && length(initial_proteins) < p2.config.run.max_initial_proteins
            p2_gene_range_size = p2_gene_range.stop - p2_gene_range.start + 1
            new_protein = deepcopy(p2.initial_cell_proteins[i])
            old_concs = new_protein.concs[p2_gene_range]
            resize!(new_protein.concs, total_concs)

            #shift the concs over the p2 genes...
            shifted_range = p1_gene_range.stop + 1 : (p1_gene_range.stop + 1) + p2_gene_range_size - 1
            new_protein.concs[shifted_range] .= old_concs
            #...and zero out the concs over the p1 genes
            new_protein.concs[1:shifted_range.start] .= 0
            new_protein.concs[shifted_range.stop + 1:total_concs] .= 0

            push!(initial_proteins, new_protein)
        end
    end

    config = Config(p1.config.run, UInt64(p1.config.rng.seed[1]) + UInt64(seed_offset))
    Individual(config, genes, initial_proteins)
end

end

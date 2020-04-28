module GeneStateMod

using GeneMod
using ProteinMod
using ProteinPropsMod
using RunMod
using MiscUtilsMod

import Base.show

export GeneState

mutable struct GeneState
    config::Config
    gene::Gene
    bindings::Array{Union{Protein, Nothing}, 1}

    function GeneState(config::Config, gene::Gene)
        new(
            config,
            gene,
            repeat([nothing], config.run.bind_sites_per_gene)
        )
    end
end

function show(io::IO, gs::GeneState, ilevel::Int64=0)
    iprintln(io, "GeneState", ilevel)

    iprintln(io, "genome_index: $(gs.gene.genome_index)", ilevel + 1)

    iprintln(io, "bindings", ilevel + 1)
    site_strs = Array{String, 1}()
    for protein in gs.bindings
        if protein == nothing
            push!(site_strs, "(nothing)")
        else
            push!(site_strs, ProteinPropsMod.to_str(protein.props))
        end
    end
    iprintln(io, join(site_strs, ", "), ilevel + 2)
    print(io, "\n")
end

function get_compatible_bind_site_indices(gs::GeneState, protein::Protein)
    col = gs.gene.genome_index
    thresh = gs.gene.threshold

    indices = Array{Int64, 1}()
    for i in 1:length(gs.gene.bind_sites)
        site = gs.gene.bind_sites[i]
        if (protein.type == site.type &&
            protein.loc == site.loc &&
            protein.concs[col] >= thresh)
            push!(indices, i)
        end
    end

    indices
end

#binding consumes some of the bound protein
#note: this must be called *after* produce()
#since it may lower the bound protein's conc below the binding threshold
function run_binding_consum(gs::GeneState)
    #binding consumes some of the protein
    for site_index in 1:length(gs.bindings)
        protein = gs.bindings[site_index]
        if protein != nothing
            rate = gs.gene.bind_sites[site_index].consum_rate
            col = gs.gene.genome_index
            protein.concs[col] = clamp(protein.concs[col] - rate, 0.0, 1.0)
        end
    end
end

function bind(gs::GeneState, protein::Protein, site_index::Int64)
    gs.bindings[site_index] = protein
end

function unbind(gs::GeneState, site_index::Int64)
    gs.bindings[site_index] = nothing
end

function get_binding(gs::GeneState, site_index::Int64)
    gs.bindings[site_index]
end

#Note: assumes proteins have already been bound
function get_prod_rates(gs::GeneState)
    rates = Array{NamedTuple{(:prod_index, :rate), Tuple{Int64, Float64}}, 1}()
    logic = gs.gene.bind_logic
    col = gs.gene.genome_index

    #each site independently activates the prod site at the same index
    #rate is determined by the amount the binding protein's conc exceeds the gene's binding threshold
    if logic == GeneMod.Id
        for i in 1:length(gs.bindings)
            protein = gs.bindings[i]
            if protein != nothing && protein.props.fcn == ProteinPropsMod.Activate
                protein_conc = protein.concs[col]
                threshold = gs.gene.bind_sites[i].threshold
                excess = protein_conc - threshold #will be >= 0, since it's already bound

                #scale excess from [0.0, 1.0 - threshold] to [0.0, run.max_prod_rate]
                rate = excess * (gs.config.run.max_prod_rate / (1 - threshold))
                push!(rates, (prod_index=i, rate=rate))
            end
        end

        #if all sites are active, first prod site is activated
        #rate is determined by the average excess across all binding sites    
    elseif logic == GeneMod.And
        sum = 0 #sum of excess
        max_sum = 0 #max possible excess
        i = 1
        while i <= length(gs.bindings) && gs.bindings[i] != nothing && gs.bindings[i].props.fcn == ProteinPropsMod.Activate
            threshold = gs.gene.bind_sites[i].threshold
            protein_conc = gs.bindings[i].concs[col]
            excess = protein_conc - threshold
            sum += excess
            max_sum += 1.0 - threshold
            i += 1
        end
        
        if i > length(gs.bindings) #if all sites were active
            avg = sum / length(gs.bindings)
            max_avg = max_sum / length(gs.bindings) #max possible avg
            #scale avg from [0.0, 1.0 - max_avg] to [0.0, run.max_prod_rate]
            rate = excess * (gs.config.run.max_prod_rate / (1 - max_avg))
            push!(rates, (prod_index=1, rate=rate))
        end

        #if at least one site is active, first prod site is activated
        #rate is determined by the average excess across all binding sites that are active    
    elseif logic == GeneMod.Or
        count = 0
        sum = 0
        max_sum = 0
        for i in 1:length(gs.bindings)
            if gs.bindings[i] != nothing && gs.bindings[i].props.fcn == ProteinPropsMod.Activate
                count += 1
                protein_conc = gs.bindings[i].concs[col]
                threshold = gs.gene.bind_sites[i].threshold
                excess = protein_conc - threshold
                sum += excess
                max_sum += 1.0 - threshold
            end
            i += 1
        end

        if count > 0
            avg = sum / count
            max_avg = max_sum / count
            #scale it from [0.0, 1.0 - max_avg] to [0.0, run.max_prod_rate]
            rate = excess * (gs.config.run.max_prod_rate / (1 - max_avg))
            push!(rates, (prod_index=1, rate=rate))
        end

        #if exactly one site is active, first prod site is activated
        #rate is determined by the excess on the binding site that is active
    elseif logic == GeneMod.Xor
        count = 0
        excess = 0
        max_excess = 0
        i = 1
        while i <= length(gs.bindings) && count < 2
            protein = gs.bindings[i]
            if protein != nothing && protein.props.fcn == ProteinPropsMod.Activate
                count += 1
                protein_conc = protein.concs[col]
                threshold = gs.gene.bind_sites[i].threshold
                excess = protein_conc - threshold
                max_excess = 1.0 - threshold
            end
            i += 1
        end

        if count == 1
            #scale excess from [0.0, 1.0 - max_excess] to [0.0, run.max_prod_rate]
            rate = excess * (gs.config.run.max_prod_rate / (1 - max_excess))
            push!(rates, (prod_index=1, rate=rate))
        end
    end

    rates
end

end

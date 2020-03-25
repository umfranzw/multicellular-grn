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
            repeat([nothing], config.run.num_bind_sites)
        )
    end
end

function show(io::IO, gs::GeneState, ilevel::Int64=0)
    iprintln(io, "GeneState", ilevel)

    iprintln(io, "genome_index: $(gs.gene.genome_index)", ilevel + 1)

    iprintln(io, "bindings", ilevel + 1)
    site_strs = Array{String, 1}()
    for protein in gs.site_bindings
        if site == nothing
            push!(site_strs, "(nothing)")
        else
            push!(ProteinPropsMod.to_str(protein.props))
        end
    end
    iprintln(io, join(site_strs, ", "), ilevel + 2)
    println(io, "")
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
    rates = Array{Union{Float64, Nothing}, 1}()
    logic = gs.gene.bind_logic
    threshold = gs.gene.threshold
    col = gs.gene.genome_index

    #each site independently activates the prod site at the same index
    #rate is determined by the amount the binding protein's conc exceeds the gene's binding threshold
    if logic == Gene.Id
        for site_index in 1:length(gs.bindings)
            rate = 0
            if protein != nothing
                protein_conc = protein.concs[col]
                excess = protein_conc - thresh #will be >= 0, since it's already bound

                #scale excess from [0.0, 1.0 - threshold] to [0.0, run.max_prod_rate]
                rate = excess * (config.run.max_prod_rate / (1 - threshold))
            end
            push!(rates, (prod_index=site_index, rate=rate))
        end

    #if all sites are active, first prod site is activated
    #rate is determined by the average excess across all binding sites    
    elseif logic == Gene.And
        sum = 0
        i = 1
        while i <= length(gs.bindings) && gs.bindings[i] != nothing
            protein_conc = protein.concs[col]
            excess = protein_conc - thresh
            sum += excess
            i += 1
        end
        
        if i > length(gs.bindings) #if all sites were active
            avg = sum / length(gs.bindings)
            #scale it from [0.0, 1.0 - threshold] to [0.0, run.max_prod_rate]
            rate = excess * (config.run.max_prod_rate / (1 - threshold))
        else
            rate = 0
        end
        push!(rates, (prod_index=1, rate=rate))

    #if at least one site is active, first prod site is activated
    #rate is determined by the average excess across all binding sites that are active    
    elseif logic == Gene.Or

    #if exactly one site is active, first prod site is activated
    #rate is determined by the excess on the binding site that is active
    elseif logic == Gene.Xor
        
    end

    rates
end

end

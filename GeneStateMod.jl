module GeneStateMod

using GeneMod
using ProteinMod
using ProteinPropsMod
using RunMod
using MiscUtilsMod

import Base.show

export GeneState

mutable struct GeneState
    run::Run
    gene::Gene
    bindings::Dict{Union{Type{BindSite}, Type{ProdSite}}, Array{Union{Protein, Nothing}, 1}}

    function GeneState(run::Run, gene::Gene)
        bindings = Dict{Union{Type{BindSite}, Type{ProdSite}}, Array{Union{Protein, Nothing}, 1}}()
        bindings[BindSite] = repeat([nothing], length(gene.bind_sites))
        bindings[ProdSite] = repeat([nothing], length(gene.prod_sites))
        
        new(
            run,
            gene,
            bindings
        )
    end
end

function show(io::IO, gs::GeneState, ilevel::Int64=0)
    iprintln(io, "GeneState", ilevel)

    iprintln(io, "genome_index: $(gs.gene.genome_index)", ilevel + 1)

    iprintln(io, "bindings", ilevel + 1)
    for type in keys(gs.bindings)
        iprintln(io, string(type), ilevel + 2)
        site_strs = Array{String, 1}()
        for protein in gs.bindings[type]
            if protein == nothing
                push!(site_strs, "(nothing)")
            else
                push!(site_strs, ProteinPropsMod.to_str(protein.props, protein.is_initial))
            end
        end
        iprintln(io, join(site_strs, ", "), ilevel + 3)
    end
    print(io, "\n")
end

# function get_compatible_bind_site_indices(gs::GeneState, protein::Protein)
#     col = gs.gene.genome_index
#     thresh = gs.gene.threshold

#     indices = Array{Int64, 1}()
#     for i in 1:length(gs.gene.bind_sites)
#         site = gs.gene.bind_sites[i]
#         if (protein.type == site.type &&
#             protein.loc == site.loc &&
#             protein.concs[col] >= thresh)
#             push!(indices, i)
#         end
#     end

#     indices
# end

#binding consumes some of the bound protein
#note: this must be called *after* produce()
#since it may lower the bound protein's conc below the binding threshold
function run_binding_consum(gs::GeneState)
    #binding consumes some of the protein
    for type in keys(gs.bindings)
        for site_index in 1:length(gs.bindings[type])
            protein = gs.bindings[type][site_index]
            if protein != nothing
                if type == BindSite
                    sites = gs.gene.bind_sites
                else
                    sites = gs.gene.prod_sites
                end
                rate = sites[site_index].consum_rate
                col = gs.gene.genome_index
                protein.concs[col] = clamp(protein.concs[col] - rate, 0.0, 1.0)
            end
        end
    end
end

function clear_all_bindings(gs::GeneState)
    for site_type in keys(gs.bindings)
        for i in 1:length(gs.bindings[site_type])
            gs.bindings[site_type][i] = nothing
        end
    end
end

function bind(gs::GeneState, protein::Protein, site_type::Union{Type{BindSite}, Type{ProdSite}}, site_index::Int64)
    gs.bindings[site_type][site_index] = protein
end

function unbind(gs::GeneState, site_type::Union{Type{BindSite}, Type{ProdSite}}, site_index::Int64)
    gs.bindings[site_type][site_index] = nothing
end

function get_binding(gs::GeneState, site_type::Union{Type{BindSite}, Type{ProdSite}}, site_index::Int64)
    gs.bindings[site_type][site_index]
end

#Note: assumes proteins have already been bound
function get_prod_rates(gs::GeneState)
    rates = Array{NamedTuple{(:prod_index, :rate), Tuple{Int64, Float64}}, 1}()
    logic = gs.gene.bind_logic
    col = gs.gene.genome_index

    #each site independently activates the prod site at the same index
    #rate is determined by the amount the binding protein's conc exceeds the gene's binding threshold
    if logic == GeneMod.Id
        for i in 1:length(gs.bindings[BindSite])
            protein = gs.bindings[BindSite][i]
            if protein != nothing && gs.bindings[ProdSite][i] == nothing
                protein_conc = protein.concs[col]
                threshold = gs.gene.bind_sites[i].threshold
                excess = protein_conc - threshold #will be >= 0, since it's already bound

                excess = clamp(excess * 2, 0, 1 - threshold)
                #scale excess from [0.0, 1.0 - threshold] to [0.0, run.max_prod_rate]
                rate = excess * (gs.run.max_prod_rate / (1 - threshold))
                push!(rates, (prod_index=i, rate=rate))
            end
        end

        #if all sites are active, first prod site is activated
        #rate is determined by the average excess across all binding sites
    elseif logic == GeneMod.And
        if gs.bindings[ProdSite][1] == nothing
            sum = 0 #sum of excess
            max_sum = 0 #max possible excess
            i = 1
            while i <= length(gs.bindings[BindSite]) && gs.bindings[BindSite][i] != nothing
                threshold = gs.gene.bind_sites[i].threshold
                protein_conc = gs.bindings[BindSite][i].concs[col]
                excess = protein_conc - threshold
                sum += excess
                max_sum += 1.0 - threshold
                i += 1
            end
            
            if i > length(gs.bindings[BindSite]) #if all sites were active
                avg = sum / length(gs.bindings[BindSite])
                max_avg = max_sum / length(gs.bindings[BindSite]) #max possible avg
                #scale avg from [0.0, 1.0 - max_avg] to [0.0, run.max_prod_rate]
                rate = excess * (gs.run.max_prod_rate / (1 - max_avg))
                push!(rates, (prod_index=1, rate=rate))
            end
        end

        #if at least one site is active, first prod site is activated
        #rate is determined by the average excess across all binding sites that are active
    elseif logic == GeneMod.Or
        if gs.bindings[ProdSite][1] == nothing
            count = 0
            sum = 0
            max_sum = 0
            for i in 1:length(gs.bindings[BindSite])
                if gs.bindings[BindSite][i] != nothing
                    count += 1
                    protein_conc = gs.bindings[BindSite][i].concs[col]
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
                rate = excess * (gs.run.max_prod_rate / (1 - max_avg))
                push!(rates, (prod_index=1, rate=rate))
            end
        end

        #if exactly one site is active, first prod site is activated
        #rate is determined by the excess on the binding site that is active
    elseif logic == GeneMod.Xor
        if gs.bindings[ProdSite][1] == nothing
            count = 0
            excess = 0
            max_excess = 0
            i = 1
            while i <= length(gs.bindings[BindSite]) && count < 2
                protein = gs.bindings[BindSite][i]
                if protein != nothing
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
                rate = excess * (gs.run.max_prod_rate / (1 - max_excess))
                push!(rates, (prod_index=1, rate=rate))
            end
        end
    end

    rates
end

end

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
using RegSimInfoMod
using FitnessInfoMod

import Random
import RandUtilsMod
import Base.show
import Random

export Individual,
    rand_init, run_bind

mutable struct Individual
    config::Config
    genes::Array{Gene, 1}
    cell_tree::CellTree
    initial_cell_proteins::Array{Protein, 1}
    #note: this is a value in [0.0, 1.0], where 0.0 is optimal
    fitness::Float64
    reg_sim_info::RegSimInfo #holds info about the *last* reg sim, if any
    fitness_info::Union{FitnessInfo, Nothing}
    id::Union{UInt64, Nothing} #for unique identification in the vis UI
    last_mod::Int64 #ea_step of last *genetic* (not regulatory) modification

    function Individual(config::Config, gene::Gene)
        genes = [gene]
        root_cell = Cell(config, genes)
        initial_proteins = Array{Protein, 1}([make_initial_protein(config, gene, root_cell)])
        
        Individual(config, genes, initial_proteins, root_cell)
    end

    function Individual(config::Config, genes::Array{Gene, 1}, initial_proteins::Array{Protein, 1})
        root_cell = Cell(config, genes)

        Individual(config, genes, initial_proteins, root_cell)
    end

    function Individual(config::Config, genes::Array{Gene, 1}, initial_proteins::Array{Protein, 1}, root_cell::Cell)
        cell_tree = CellTree(root_cell)
        reg_sim_info = RegSimInfo(1)

        indiv = new(config, genes, cell_tree, initial_proteins, 1.0, reg_sim_info, nothing, nothing, -1)
        indiv.id = hash(indiv)
        CellMod.insert_initial_proteins(root_cell, indiv.initial_cell_proteins)

        indiv
    end
end

function rand_init(run::Run, seed_offset::UInt64)
    config = Config(run, seed_offset)

    #genes, initial_protein_props = make_initial_genes(config)
    #genes = map(i -> GeneMod.rand_init(config, i, [ProteinPropsMod.Internal], [GeneMod.Id]), 1:config.run.num_initial_genes)

    gene = make_initial_gene(config)

    Individual(config, gene)
end

function get_id_str(indiv::Individual)
    "$(indiv.id):$(indiv.last_mod)"
end

function reset_reg_sim_info(indiv::Individual)
    RegSimInfoMod.reset(indiv.reg_sim_info)
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
        threshold=[config.run.bind_threshold],
        consum_rate=[config.run.bind_consum_rate]
    )

    ab_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Internal],
        tag=[UInt8(1)],
        fcn=ProteinPropsMod.Activate,
        threshold=[config.run.bind_threshold],
        consum_rate=[config.run.bind_consum_rate]
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
        threshold=[config.run.bind_threshold],
        consum_rate=[config.run.bind_consum_rate]
    )

    bc_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Application],
        tag=[UInt8(2)],
        action=[ProteinPropsMod.Divide],
        threshold=[config.run.bind_threshold],
        consum_rate=[config.run.bind_consum_rate]
    )
    bc = pop_remaining_sites(config, genome_index, bc_bind_site, bc_prod_site)

    #ba
    ba_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ab_prod_site.type],
        tag=[ab_prod_site.tag],
        threshold=[config.run.bind_threshold],
        consum_rate=[config.run.bind_consum_rate]
    )

    ba_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ab_bind_site.type],
        tag=[ab_bind_site.tag],
        fcn=ProteinPropsMod.Activate,
        threshold=[config.run.bind_threshold],
        consum_rate=[config.run.bind_consum_rate]
    )
    ba = pop_remaining_sites(config, genome_index, ba_bind_site, ba_prod_site)

    #de
    de_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ProteinPropsMod.Internal],
        tag=[UInt8(3)],
        threshold=[config.run.bind_threshold],
        consum_rate=[config.run.bind_consum_rate]
    )

    plus_index = findfirst(i -> SymProbsMod.index_to_sym[i].val == :+, 1:length(SymProbsMod.index_to_sym)) - 1
    de_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Application],
        tag=[UInt8(4)],
        action=[ProteinPropsMod.SymProb],
        fcn=ProteinPropsMod.Activate,
        arg=[Int8(plus_index)],
        threshold=[config.run.bind_threshold],
        consum_rate=[config.run.bind_consum_rate]
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
        threshold=[config.run.bind_threshold],
        consum_rate=[config.run.bind_consum_rate]
    )

    bf_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Neighbour],
        tag=[UInt8(4)],
        fcn=ProteinPropsMod.Activate,
        arg=[Int8(4)], #right child
        threshold=[config.run.bind_threshold],
        consum_rate=[config.run.bind_consum_rate]
    )
    bf = pop_remaining_sites(config, genome_index, bf_bind_site, bf_prod_site)

    #bg
    bg_bind_site = GeneMod.rand_bind_site(
        config,
        type=[ab_prod_site.type],
        tag=[ab_prod_site.tag],
        threshold=[config.run.bind_threshold],
        consum_rate=[config.run.bind_consum_rate]
    )

    bg_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Neighbour],
        tag=[UInt8(5)],
        fcn=ProteinPropsMod.Activate,
        arg=[Int8(3)], #left child
        threshold=[config.run.bind_threshold],
        consum_rate=[config.run.bind_consum_rate]
    )
    bg = pop_remaining_sites(config, genome_index, bg_bind_site, bg_prod_site)

    #gh
    gh_bind_site = GeneMod.rand_bind_site(
        config,
        type=[bg_prod_site.type],
        tag=[bg_prod_site.tag],
        threshold=[config.run.bind_threshold],
        consum_rate=[config.run.bind_consum_rate]
    )

    x_index = findfirst(i -> SymProbsMod.index_to_sym[i].val == :x, 1:length(SymProbsMod.index_to_sym)) - 1
    gh_prod_site = GeneMod.rand_prod_site(
        config,
        type=[ProteinPropsMod.Application],
        tag=[UInt8(6)],
        action=[ProteinPropsMod.SymProb],
        fcn=ProteinPropsMod.Activate,
        arg=[Int8(x_index)],
        threshold=[config.run.bind_threshold],
        consum_rate=[config.run.bind_consum_rate]
    )
    gh = pop_remaining_sites(config, genome_index, gh_bind_site, gh_prod_site)

#fnc
fnc_bind_site = GeneMod.rand_bind_site(
    config,
    type=[bf_prod_site.type], #neighbour-only receptor
    tag=[bf_prod_site.tag],
    threshold=[config.run.bind_threshold],
    consum_rate=[config.run.bind_consum_rate]
)

fnc_prod_site = GeneMod.rand_prod_site(
    config,
    type=[ProteinPropsMod.Neighbour],
    tag=[bc_prod_site.tag],
    fcn=ProteinPropsMod.Inhibit,
    arg=[Int8(1)], #parent
    threshold=[config.run.bind_threshold],
    consum_rate=[config.run.bind_consum_rate]
)
fnc = pop_remaining_sites(config, genome_index, fnc_bind_site, fnc_prod_site)

#fi
fi_bind_site = GeneMod.rand_bind_site(
    config,
    type=[bf_prod_site.type], #neighbour-only receptor
    tag=[bf_prod_site.tag],
    threshold=[config.run.bind_threshold],
    consum_rate=[config.run.bind_consum_rate]
)

one_index = findfirst(i -> SymProbsMod.index_to_sym[i].val == :1, 1:length(SymProbsMod.index_to_sym)) - 1
fi_prod_site = GeneMod.rand_prod_site(
    config,
    type=[ProteinPropsMod.Application],
    tag=[UInt8(7)],
    action=[ProteinPropsMod.SymProb],
    fcn=ProteinPropsMod.Activate,
    arg=[Int8(one_index)],
    threshold=[config.run.bind_threshold],
    consum_rate=[config.run.bind_consum_rate]
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
        push!(bind_sites, GeneMod.rand_bind_site(config, consum_rate=[config.run.bind_consum_rate],  threshold=[config.run.bind_threshold]))
        push!(prod_sites, GeneMod.rand_prod_site(config, consum_rate=[config.run.bind_consum_rate],  threshold=[config.run.bind_threshold]))
    end

    Gene(config, genome_index, GeneMod.Id, bind_sites, prod_sites)
end

function make_initial_gene(config::Config)
    GeneMod.rand_init(
        config,
        1,
        bind_site_types=[ProteinPropsMod.Internal],
        bind_logic=[GeneMod.Id]
    )
end

function make_initial_protein(config::Config, initial_gene::Gene, root_cell::Cell)
    half = (1.0 - config.run.bind_threshold) / 2

    props = ProteinPropsMod.rand_init(
        config,
        type=[ProteinPropsMod.Internal],
        tag=[initial_gene.bind_sites[1].tag],
        fcn=ProteinPropsMod.Activate
    )
    
    protein = Protein(config, props, false, true, 1, root_cell.id)
    protein.concs = RandUtilsMod.rand_floats(config, config.run.bind_threshold + half, 1.0, 1)
    
    protein
end

#resets everything to the way it was before the reg sim (so
#we can run the reg sim again on the next ea_step)
function reset_cell_tree(indiv::Individual)
    #just re-initialize the cell (this discards the rest of the tree, along with any protein bindings)
    indiv.cell_tree.root = Cell(indiv.config, indiv.genes)
    CellMod.insert_initial_proteins(indiv.cell_tree.root, indiv.initial_cell_proteins)
end

function show(io::IO, indiv::Individual, ilevel::Int64=0)
    iprintln(io, "Individual:", ilevel)
    iprintln(io, "genes:", ilevel + 1)
    foreach(g -> GeneMod.show(io, g, ilevel + 2), indiv.genes)

    iprintln(io, "cell_tree:", ilevel + 1)
    iprintln(io, indiv.cell_tree, ilevel + 2)

    iprintln(io, "initial_cell_proteins:", ilevel + 1)
    foreach(p -> ProteinMod.show(io, p, ilevel + 2), indiv.initial_cell_proteins)

    iprintln(io, "fitness: $(indiv.fitness)", ilevel + 1)
end

function run_fix_syms(indiv::Individual)
    CellTreeMod.traverse(CellMod.fix_sym, indiv.cell_tree)
end

function count_proteins(indiv::Individual)
    num_proteins = 0
    CellTreeMod.traverse(c -> num_proteins += ProteinStoreMod.num_proteins(c.proteins), indiv.cell_tree)

    num_proteins
end

function has_proteins(indiv::Individual)
    has_proteins(indiv.cell_tree.root)
end

function has_proteins(cell::Cell)
    result = ProteinStoreMod.has_proteins(cell.proteins)
    
    i = 1
    while !result && i <= length(cell.children)
        result = has_proteins(cell.children[i])
        i += 1
    end

    result
end

end

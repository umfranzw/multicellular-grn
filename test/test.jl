using RunMod
using IndividualMod
using GeneMod
using GeneStateMod
using ProteinPropsMod
using CellMod
using CellTreeMod
using SymProbsMod
using ProteinStoreMod
using ProteinMod
using RandUtilsMod

import Random

function get_config()
    run = RunMod.get_first_run()
    if run.fix_rng_seed
        seed_base = run.rng_seed
    else
        dev = Random.RandomDevice()
        seed_base = UInt64(Random.rand(dev) * 0xffffffffffffffff)
    end
    config = Config(run, Random.MersenneTwister(seed_base))
end

function get_i1(config::Config)
    bind_sites = Array{BindSite, 1}()
    prod_sites = Array{ProteinProps, 1}()
    for i in 1:config.run.bind_sites_per_gene
        push!(bind_sites, BindSite(
            ProteinPropsMod.Internal,
            ProteinPropsMod.SymProb,
            ProteinPropsMod.Bottom,
            0.1,
            0.1
        ))
        push!(prod_sites, ProteinProps(
            ProteinPropsMod.Application,
            ProteinPropsMod.Activate,
            ProteinPropsMod.Sensor,
            ProteinPropsMod.Bottom,
            UInt8(10 + i)
        ))
    end

    gene = Gene(
        config,
        1,
        GeneMod.Id,
        bind_sites,
        prod_sites
    )
    genes = [gene]

    sensors = Dict{ProteinPropsMod.ProteinLoc, Array{Float64, 1}}()
    for loc in instances(ProteinPropsMod.ProteinLoc)
        sensors[loc] = [config.run.initial_cell_sensor_conc]
    end

    root = Cell(
        config,
        genes
    )

    initial_proteins = Array{Protein, 1}()
    for i in 1:config.run.num_initial_proteins
        gene = genes[(i % length(genes)) + 1]
        types = Array{ProteinPropsMod.ProteinType, 1}()
        actions = Array{ProteinPropsMod.ProteinAction, 1}()
        locs = Array{ProteinPropsMod.ProteinLoc, 1}()
        max_thresh = 0.0
        for site in gene.bind_sites
            push!(types, site.type)
            push!(actions, site.action)
            push!(locs, site.loc)
            max_thresh = max(site.threshold, max_thresh)
        end
        
        props = ProteinPropsMod.rand_init(
            config,
            type=types,
            fcn=[ProteinPropsMod.Activate],
            action=actions,
            loc=locs
        )

        #note: it is possible that not all initial proteins in the array are unique. That's ok, since they'll be subject to evolution.
        protein = Protein(config, props, false, true, length(root.gene_states), pointer_from_objref(root))
        #protein.concs = RandUtilsMod.rand_floats(config, max_thresh, 1.0, length(genes))
        protein.concs = repeat([1.0], length(genes))
        push!(initial_proteins, protein)
    end

    indiv = Individual(
        config,
        [gene],
        CellTree(root),
        initial_proteins,
        1.0,
        Array{Int64, 1}([0])
    )
    CellMod.insert_initial_proteins(root, initial_proteins)

    indiv
end

function get_i1(config::Config)
    bind_sites1 = Array{BindSite, 1}()
    prod_sites1 = Array{ProteinProps, 1}()
    for i in 1:config.run.bind_sites_per_gene
        push!(bind_sites1, BindSite(
            ProteinPropsMod.Internal,
            ProteinPropsMod.SymProb,
            ProteinPropsMod.Bottom,
            0.1,
            0.1
        ))
        push!(prod_sites1, ProteinProps(
            ProteinPropsMod.Neighbour,
            ProteinPropsMod.Activate,
            ProteinPropsMod.Sensor,
            ProteinPropsMod.Bottom,
            UInt8(10 + i)
        ))
    end

    gene1 = Gene(
        config,
        1,
        GeneMod.Id,
        bind_sites1,
        prod_sites1
    )

    bind_sites2 = Array{BindSite, 1}()
    prod_sites2 = Array{ProteinProps, 1}()
    for i in 1:config.run.bind_sites_per_gene
        push!(bind_sites2, BindSite(
            ProteinPropsMod.Neighbour,
            ProteinPropsMod.Sensor,
            ProteinPropsMod.Top,
            0.1,
            0.1
        ))
        push!(prod_sites2, ProteinProps(
            ProteinPropsMod.Application,
            ProteinPropsMod.Activate,
            ProteinPropsMod.SymProb,
            ProteinPropsMod.Bottom,
            UInt8(10 + i)
        ))
    end
    
    gene2 = Gene(
        config,
        1,
        GeneMod.Id,
        bind_sites2,
        prod_sites2
    )
    
    genes = [gene1, gene2]

    sensors = Dict{ProteinPropsMod.ProteinLoc, Array{Float64, 1}}()
    for loc in instances(ProteinPropsMod.ProteinLoc)
        sensors[loc] = [config.run.initial_cell_sensor_conc]
    end

    root = Cell(
        config,
        genes
    )

    initial_proteins = Array{Protein, 1}()
    for i in 0:config.run.num_initial_proteins - 1
        gene = genes[(i % length(genes)) + 1]
        types = Array{ProteinPropsMod.ProteinType, 1}()
        actions = Array{ProteinPropsMod.ProteinAction, 1}()
        locs = Array{ProteinPropsMod.ProteinLoc, 1}()
        max_thresh = 0.0
        for site in gene.bind_sites
            push!(types, site.type)
            push!(actions, site.action)
            push!(locs, site.loc)
            max_thresh = max(site.threshold, max_thresh)
        end
        
        props = ProteinPropsMod.rand_init(
            config,
            type=types,
            fcn=[ProteinPropsMod.Activate],
            action=actions,
            loc=locs
        )

        #note: it is possible that not all initial proteins in the array are unique. That's ok, since they'll be subject to evolution.
        protein = Protein(config, props, false, true, length(root.gene_states), pointer_from_objref(root))
        #protein.concs = RandUtilsMod.rand_floats(config, max_thresh, 1.0, length(genes))
        protein.concs = repeat([1.0], length(genes))
        push!(initial_proteins, protein)
    end

    indiv = Individual(
        config,
        genes,
        CellTree(root),
        initial_proteins,
        1.0,
        repeat([0], length(genes))
    )
    CellMod.insert_initial_proteins(root, initial_proteins)

    child = Cell(
        config,
        genes
    )
    CellMod.add_parent(child, root)
    CellMod.insert_initial_proteins(child, initial_proteins)

    indiv
end

config = get_config()
i1 = get_i1(config)
root = i1.cell_tree.root
child = root.children[1]

IndividualMod.run_bind(i1)
IndividualMod.run_produce(i1)
IndividualMod.run_binding_consum(i1)
IndividualMod.run_diffuse(i1)
# IndividualMod.run_neighbour_comm(indiv)
# childps = map(identity, ProteinStoreMod.get_by_type(child.proteins, ProteinPropsMod.Neighbour))
# rootps = map(identity, ProteinStoreMod.get_by_type(root.proteins, ProteinPropsMod.Neighbour))

module GeneMod

using RunMod
using ProteinMod

import RandUtilsMod
import Random

export Gene,
    rand_init

mutable struct Gene
    run::Run
    genome_index::Int64
    reg_site::BitArray{1}
    growth_site::BitArray{1}
    bind_sites::Array{BitArray{1}, 1}
    prod_sites::Array{BitArray{1}, 1}
end

function extend_rand(run::Run, bits::BitArray{1}, n)
    if length(bits) < n
        needed = n - length(bits)
        bits = cat(bits, RandUtilsMod.rand_bits(run, needed), dims=1)
    end

    bits
end

function rand_reg_site(run::Run, genome_index::Int64)
    scope = Random.rand(instances(ProteinMod.Scope))
    target = ProteinMod.Internal
    affinity = ProteinMod.Reg
    action = Random.rand(instances(ProteinMod.RegAction))

    bits = BitArray(Int64(scope), min_bits=ProteinMod.num_scope_bits)
    bits = cat(bits, BitArray(Int64(target), min_bits=ProteinMod.num_target_bits), dims=1)
    bits = cat(bits, BitArray(Int64(affinity), min_bits=ProteinMod.num_affinity_bits), dims=1)
    bits = cat(bits, BitArray(Int64(action), min_bits=ProteinMod.num_reg_action_bits), dims=1)
    
    extend_rand(run, bits, ProteinMod.num_bits)
end

function rand_growth_site(run::Run, genome_index::Int64)
    scope = Random.rand(instances(ProteinMod.Scope))
    target = ProteinMod.Internal
    affinity = ProteinMod.Grow
    dir = Random.rand(instances(ProteinMod.GrowthDir))

    bits = BitArray(Int64(scope), min_bits=ProteinMod.num_scope_bits)
    bits = cat(bits, BitArray(Int64(target), min_bits=ProteinMod.num_target_bits), dims=1)
    bits = cat(bits, BitArray(Int64(affinity), min_bits=ProteinMod.num_affinity_bits), dims=1)
    bits = cat(bits, BitArray(Int64(dir), min_bits=ProteinMod.num_growth_dir_bits), dims=1)

    extend_rand(run, bits, ProteinMod.num_bits)
end

function rand_bind_site(run::Run, genome_index::Int64)
    scope = Random.rand(instances(ProteinMod.Scope))
    target = ProteinMod.Internal
    affinity = Random.rand([ProteinMod.Bind])

    bits = BitArray([Bool(scope)])
    bits = cat(bits, BitArray([Bool(target)]), dims=1)
    bits = cat(bits, BitArray(Int64(affinity)), dims=1)
    
    extend_rand(run, bits, ProteinMod.num_bits)
end

function rand_prod_site(run::Run, genome_index::Int64)
    RandUtilsMod.rand_bits(run, ProteinMod.num_bits) #no constraints
end

function rand_init(run::Run, genome_index::Int64)
    reg_site = rand_reg_site(run, genome_index)
    growth_site = rand_growth_site(run, genome_index)

    bind_sites = []
    prod_sites = []

    for i in 1:run.num_bind_sites
        push!(bind_sites, rand_bind_site(run, genome_index))
        push!(prod_sites, rand_prod_site(run, genome_index))
    end
    
    Gene(
        run,
        genome_index,
        reg_site,
        growth_site,
        bind_sites,
        prod_sites
    )
end

end

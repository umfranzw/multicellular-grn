module ProteinMod

using RunMod
using BitUtilsMod
import RandUtilsMod
import Random
import Base.copy

export Scope, Target,
    Protein,
    get_scope, get_fcn, rand_init,
    num_bits, num_affinity_bits, num_fcn_bits

@enum Scope::Bool IntraCell=false InterCell=true
@enum Target::Bool Internal=false Output=true
#note: for now, length must be a power of 2
@enum BindAffinity::Int64 Bind=0 Prod=1 Reg=2 Grow=3
@enum RegAction::Bool RateUp=false RateDown=true
@enum GrowthDir::Bool GrowDown=false GrowUp=true

struct ProteinFcn
    name::String
    lisp_code::String
end

#note: for now, length must be a power of 2
fcns = [
    ProteinFcn("a", ""),
    ProteinFcn("b", ""),
    ProteinFcn("c", ""),
    ProteinFcn("d", "")
]

const num_scope_bits = BitUtilsMod.bits_required(length(instances(Scope)))
const num_target_bits = BitUtilsMod.bits_required(length(instances(Target)))
const num_fcn_bits = BitUtilsMod.bits_required(length(fcns))
const num_affinity_bits = BitUtilsMod.bits_required(length(instances(BindAffinity)))
const num_reg_action_bits = BitUtilsMod.bits_required(length(instances(RegAction)))
const num_growth_dir_bits = BitUtilsMod.bits_required(length(instances(GrowthDir)))

const num_bits = num_scope_bits + num_target_bits + max(num_fcn_bits, num_affinity_bits + num_reg_action_bits, num_affinity_bits + num_growth_dir_bits)

const scope_range = 1:num_scope_bits
next = scope_range.stop + 1
const target_range = next:next + num_target_bits
next = target_range.stop + 1
const fcn_range = next:next + num_fcn_bits
const affinity_range = next:next + num_affinity_bits
next = affinity_range.stop + 1
const reg_action_range = next:next + num_reg_action_bits
const growth_dir_range = next:next + num_growth_dir_bits

fcn_dict = Dict{BitArray{1}, ProteinFcn}()
for i in 0:length(fcns) - 1
    bits = BitArray(i, min_bits=num_bits)
    fcn_dict[bits] = fcns[i + 1]
end

mutable struct Protein
    run::Run
    seq::BitArray{1}
    concs::Array{Array{Float64, 1}, 1}

    function Protein(run::Run, seq::BitArray{1}, num_cells::Int64, rand_concs::Bool)
        if rand_concs
            concs = map(i -> RandUtilsMod.rand_floats(run, run.num_genes), 1:num_cells)
        else
            concs = map(i -> zeros(Float64, run.num_genes), 1:num_cells)
        end
        
        new(
            run,
            seq,
            concs
        )
    end

    function Protein(run::Run, seq::BitArray{1}, concs::Array{Array{Float64, 1}, 1})
        new(run, seq, concs)
    end
end

function copy(protein::Protein)
    #don't need to copy the run, just the other props
    Protein(protein.run, copy(protein.seq), copy(protein.concs))
end

function get_scope(protein::Protein)
    get_scope(protein.seq)
end

function get_scope(seq::BitArray{1})
    Scope(Int64(seq[scope_range]))
end

function get_target(protein::Protein)
    get_target(protein.seq)
end

function get_target(seq::BitArray{1})
    Target(Int64(seq[target_range]))
end

function get_bind_affinity(protein::Protein)
    get_bind_affinity(protein.seq)
end

function get_bind_affinity(seq::BitArray{1})
    BindAffinity(Int64(seq[affinity_range]))
end

function get_reg_action(protein::Protein)
    get_reg_action(protein.seq)
end

function get_reg_action(seq::BitArray{1})
    RegAction(Int64(seq[reg_action_range]))
end

function get_growth_dir(protein::Protein)
    get_growth_dir(protein.seq)
end

function get_growth_dir(seq::BitArray{1})
    GrowthDir(Int64(seq[growth_dir_range]))
end

function get_fcn(protein::Protein)
    get_fcn(protein.seq)
end

function get_fcn(seq::BitArray{1})
    fcn_dict[Int64(seq[fcn_range]) + 1]
end

end

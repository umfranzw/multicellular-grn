module ProteinMod

using RunMod
using BitUtilsMod
import RandUtilsMod
import Random

export Scope, Target,
    Protein,
    get_scope, get_fcn, rand_init,
    num_bits, num_affinity_bits, num_fcn_bits

@enum Scope::Bool IntraCell=false InterCell=true
@enum Target::Bool Internal=false Output=true
#note: for now, length must be a power of 2
@enum BindAffinity::Int64 Bind=0 Prod=1 RegUp=2 RegDown=3

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

const num_fcn_bits = Int64(ceil(log2(length(fcns))))
const num_affinity_bits = Int64(ceil(log2(length(instances(BindAffinity)))))
const num_info_bits = 2 #scope, target
const num_bits = num_info_bits + max(num_fcn_bits, num_affinity_bits)

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
end

function get_scope(protein::Protein)
    get_scope(protein.seq)
end

function get_scope(seq::BitArray{1})
    Scope(seq[1])
end

function get_target(protein::Protein)
    get_target(protein.seq)
end

function get_target(seq::BitArray{1})
    Target(seq[2])
end

function get_bind_affinity(protein::Protein)
    get_bind_affinity(protein.seq)
end

function get_bind_affinity(seq::BitArray{1})
    BindAffinity(Int64(seq[3:3 + num_affinity_bits - 1]))
end

function get_fcn(protein::Protein)
    get_fcn(protein.seq)
end

function get_fcn(seq::BitArray{1})
    fcn_dict[seq[3:3 + num_fcn_bits - 1]]
end

end

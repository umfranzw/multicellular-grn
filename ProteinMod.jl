module ProteinMod

using RunMod
using BitUtilsMod
import RandUtilsMod

export Scope, Effect, Target,
    Protein,
    get_scope, get_effect, get_fcn, rand_init,
    num_bits

@enum Scope::Bool IntraCell=false InterCell=true
@enum Effect::Bool Inhibit=false Activate=true
@enum Target::Bool Internal=false Output=true

struct ProteinFcn
    name::String
    lisp_code::String
end

fcns = [
    ProteinFcn("a", ""),
    ProteinFcn("b", ""),
    ProteinFcn("c", "")
]
const num_fcn_bits = Int64(ceil(log2(length(fcns))))
const num_info_bits = 3 #scope, effect, target
const num_bits = num_fcn_bits + num_info_bits

fcn_dict = Dict{BitArray{1}, ProteinFcn}()

for i in 0:length(fcns) - 1
    bits = BitArray(i, min_bits=num_bits)
    fcn_dict[bits] = fcns[i + 1]
end
    
mutable struct Protein
    run::Run
    seq::BitArray{1}
    concs::Array{Array{Float64, 1}, 1}
end

#note: does not randomly initialize concs (only seq)
function rand_init(run::Run, num_cells::Int64)
    Protein(
        run,
        RandUtilsMod.rand_bits(run, num_bits),
        map(i -> zeros(Float64, run.num_genes), 1:num_cells)
    )
end

function get_scope(protein::Protein)
    Scope(protein.seq[1])
end

function get_effect(protein::Protein)
    Effect(protein.seq[2])
end

function get_fcn(protein::Protein)
    fcn_dict[protein.seq[3:end]]
end

end

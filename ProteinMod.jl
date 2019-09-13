module ProteinMod

using RunMod
using BitUtilsMod

export ProteinScope, ProteinEffect, ProteinTarget,
    Protein,
    get_scope, get_effect, get_fcn

@enum ProteinScope::Bool IntraCell=false InterCell=true
@enum ProteinEffect::Bool Inhibit=false Activate=true
@enum ProteinTarget::Bool Internal=false Output=true

struct ProteinFcn
    name::String
    lisp_code::String
end

ProteinFcns = [
    ProteinFcn("a", "a"),
    ProteinFcn("b", "b"),
    ProteinFcn("c", "c")
]
fcn_bits = Int64(ceil(log2(length(ProteinFcns))))
ProteinFcnDict = Dict{BitArray{1}, ProteinFcn}()

for i in 0:length(ProteinFcns) - 1
    bits = BitArray(i, min_bits=fcn_bits)
    ProteinFcnDict[bits] = ProteinFcns[i + 1]
end
    
mutable struct Protein
    run::Run
    seq::BitArray{1}
    concs::Array{Float64, 1}
end

function get_scope(protein::Protein)
    ProteinScope(protein.seq[1])
end

function get_effect(protein::Protein)
    ProteinEffect(protein.seq[2])
end

function get_fcn(protein::Protein)
    ProteinFcnDict[protein.seq[3:3 + fcn_bits - 1]]
end

end

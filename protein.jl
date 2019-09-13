include("run.jl")

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
ProteinFcnDict = Dict{BitArray{1}, ProteinFcn}(zip(map(BitArray, 0:length(ProteinFcns) - 1), ProteinFcns))

struct Protein
    Run::Run
    seq::BitArray{1}
end

function get_scope(protein::Protein)
    ProteinScope(protein.seq[1])
end

function get_effect(protein::Protein)
    ProteinEffect(protein.seq[2])
end

function get_fcn(protein::Protein)
    ProteinFcnDict[protein.seq[3:7]]
end


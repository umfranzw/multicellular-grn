module ProteinMod

import RandUtilsMod
import Base.copy
import Base.show
import Formatting

using RunMod
using ProteinPropsMod
using MiscUtilsMod

export Protein,
    copy

const conc_fs = Formatting.FormatSpec("0.3f")

mutable struct Protein
    run::Run
    props::ProteinProps
    concs::Array{Float64, 1}

    function Protein(run::Run, props::ProteinProps, rand_concs::Bool)
        if rand_concs
            concs = RandUtilsMod.rand_floats(run, run.num_genes)
        else
            concs = zeros(Float64, run.num_genes)
        end
        
        new(run, props, concs)
    end

    function Protein(run::Run,  props::ProteinProps, concs::Array{Float64, 1})
        new(run, props, concs)
    end
end

function copy(protein::Protein)
    #only the concs need to be deep copied
    Protein(protein.run, protein.props, copy(protein.concs))
end

function show(io::IO, protein::Protein, ilevel::Int64=0)
    iprintln(io, "Protein:", ilevel)
    iprint(io, "props: $(protein.props)", ilevel + 1)
    
    iprint(io, "concs: [", ilevel + 1)
    for i in 1:length(protein.concs)
        Formatting.printfmt(io, conc_fs, protein.concs[i])
        if i < length(protein.concs)
            print(io, ", ")
        end
    end
    print(io, "]")
    println(io, "")
end

end

module ProteinMod

import RandUtilsMod
import Base.show
import Formatting

using RunMod
using ProteinPropsMod

export Protein

const conc_fs = Formatting.FormatSpec("0.3f")

mutable struct Protein
    props::ProteinProps
    concs::Array{Float64, 1}
    is_initial::Bool

    function Protein(config::Config, props::ProteinProps, rand_concs::Bool, is_initial::Bool, num_concs::Int64)
        if rand_concs
            concs = RandUtilsMod.rand_floats(config, num_concs)
        else
            concs = zeros(Float64, num_concs)
        end
        
        new(props, concs, is_initial)
    end
end

function show(io::IO, protein::Protein, ilevel::Int64=0)
    iprintln(io, "Protein:", ilevel)
    iprint(io, "is_initial: $(string(protein.is_initial))", ilevel + 1)
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

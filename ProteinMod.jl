module ProteinMod

import RandUtilsMod
import Base.show

using RunMod
using ProteinPropsMod
using Printf
using MiscUtilsMod

export Protein

mutable struct Protein
    props::ProteinProps
    concs::Array{Float64, 1}
    is_initial::Bool
    src_cell_id::UInt64

    function Protein(config::Config, props::ProteinProps, rand_concs::Bool, is_initial::Bool, num_concs::Int64, src_cell_id::UInt64)
        if rand_concs
            concs = RandUtilsMod.rand_floats(config, num_concs)
        else
            concs = zeros(Float64, num_concs)
        end
        
        new(props, concs, is_initial, src_cell_id)
    end
end

function show(io::IO, protein::Protein, ilevel::Int64=0)
    iprintln(io, "Protein:", ilevel)
    iprint(io, "is_initial: $(string(protein.is_initial))", ilevel + 1)
    iprint(io, "props: $(protein.props)", ilevel + 1)
    
    str_concs = join(map(c -> @sprintf("%0.2f", c), protein.concs), ", ")
    iprintln(io, "concs: [$(str_concs)]", ilevel + 1)
    iprintln(io, "src_cell_id: $(protein.src_cell_id)", ilevel + 1)
end

end

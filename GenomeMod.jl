module GenomeMod

using RunMod
using GeneMod

mutable struct Genome
    run::Run
    genes::Array{Gene, 1}
    initial_proteins::Array{Protein, 1}
    proteins::Dict{BitArray, Protein}
end

end

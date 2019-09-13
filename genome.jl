struct Genome
    run::Run
    genes::Array{Gene, 1}
    initial_proteins::Array{Protein, 1}
    proteins::Dict{BitArray, Protein}
end

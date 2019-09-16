module IndividualMod

struct Individual
    genome::Genome
    cells::Array{Cell, 1}
    initial_proteins::Array{Protein, 1}
end

end

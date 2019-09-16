module GenomeMod

using RunMod
using GeneMod

export Genome,
    rand_init

mutable struct Genome
    run::Run
    genes::Array{Gene, 1}
end

function rand_init(run::Run)
    genes = Array{Gene, 1}()
    for i in 1:run.num_genes
        push!(genes, GeneMod.rand_init(run, i))
    end
    
    Genome(run, genes)
end

end

struct Gene
    run::Run
    genome_index::Int64
    initial_output_rate::Float64
    output_rate_threshold::Float64
    binding_threshold::Float64
    binding_seqs::Array{BitArray{1}, 1}
    exons::Array{Exon{1}, 1}
end

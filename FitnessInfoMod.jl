module FitnessInfoMod

export FitnessInfo

mutable struct FitnessInfo
    contains_x::Float64
    contains_fncall::Float64
    bind_coverage::Float64
    prod_coverage::Float64
    divided::Float64
    altered_sym_prob::Float64
    genome_len::Float64
    accuracy::Float64
end

end

module FitnessInfoMod

export FitnessInfo

mutable struct FitnessInfo
    contains_x::Float64
    bind_coverage::Float64
    prod_coverage::Float64
    divided::Float64
    accuracy::Float64
end

end
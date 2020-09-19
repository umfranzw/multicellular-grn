module BestInfoMod

using IndividualMod

export BestInfo

mutable struct BestInfo
    index::Union{Nothing, Tuple{Int64, Int64, Int64}} #(ea_step, pop_index, reg_step) - reg step will always be run.reg_steps + 1
    indiv::Union{Nothing, Individual}

    function BestInfo()
        new(nothing, nothing)
    end
end

function update(info::BestInfo, indiv::Individual, ea_step::Int64, pop_index::Int64, reg_step::Int64; make_copy::Bool=true)
    updated = false
    if !is_set(info) || indiv.fitness < info.indiv.fitness
        if make_copy
            info.indiv = deepcopy(indiv)
        else
            info.indiv = indiv
        end
        
        info.index = (ea_step, pop_index, reg_step)
        updated = true
    end

    updated
end

function is_set(info::BestInfo)
    info.indiv != nothing
end

end

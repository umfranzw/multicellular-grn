module RegSimInfoMod

export RegSimInfo

mutable struct RegSimInfo
    produce_count::Array{Int64, 1}
    bind_count::Array{Int64, 1}
    division_count::Int64
    alter_sym_prob_count::Int64

    function RegSimInfo(num_genes::Int64)
        new(zeros(num_genes), zeros(num_genes), 0, 0)
    end
end

function reset(info::RegSimInfo)
    info.produce_count .= 0
    info.bind_count .= 0
    info.division_count = 0
    info.alter_sym_prob_count = 0
end

get_bind_coverage(info::RegSimInfo) = get_coverage(info.bind_count)
get_prod_coverage(info::RegSimInfo) = get_coverage(info.produce_count)
get_coverage(count::Array{Int64, 1}) = sum(count .> 0) / length(count)

end

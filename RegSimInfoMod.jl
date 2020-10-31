module RegSimInfoMod

export RegSimInfo

mutable struct RegSimInfo
    produce_count::Array{Int64, 1}
    bind_count::Array{Int64, 1}
    division_count::Int64
    alter_sym_prob_count::Int64
    reg_step_count::Int64

    function RegSimInfo(num_genes::Int64)
        new(zeros(num_genes), zeros(num_genes), 0, 0, 0)
    end
end

function reset(info::RegSimInfo)
    info.produce_count .= 0
    info.bind_count .= 0
    info.division_count = 0
    info.alter_sym_prob_count = 0
    info.reg_step_count = 0
end

function insert_new_counts(info::RegSimInfo, index::Int64, val::Int64=0)
    insert!(info.produce_count, index, val)
    insert!(info.bind_count, index, val)
end

get_bind_coverage(info::RegSimInfo) = get_coverage(info.bind_count)
get_prod_coverage(info::RegSimInfo) = get_coverage(info.produce_count)
get_coverage(count::Array{Int64, 1}) = sum(count .> 0) / length(count)
divided(info::RegSimInfo) = info.division_count > 0
altered_sym_probs(info::RegSimInfo) = info.alter_sym_prob_count > 0
get_avg_binds_per_gene(info::RegSimInfo) = info.bind_count / length(info.bind_count)

end

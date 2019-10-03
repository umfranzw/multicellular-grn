module MiscUtilsMod

function num_enum_vals(enum::Any)
    length(instances(enum))
end

function digits_needed(n::Int64)
    Int64(ceil(log10(n)))
end

end

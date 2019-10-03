module ExprUtils

import Espresso

function find_all_vars(ex::Expr, names::Array{Symbol, 1})
    vars = Array{Union{SymRef, ExprRef}, 1}()

    for name in names
        result = Espresso.findex(name, ex)

        for i in 1:length(result)
            cur = result[i]
            if cur isa Expr
                push!(vars, cur)
                for arg in cur.args[2:end]
                    push!(items, find_all_vars(arg, names)...)
                end
            else
                push!(items, ex)
            end
        end
    end

    items
end

function find_all_refs(ex::Expr)
    items = []
    pat = :(_a[_i])

    result = Espresso.findex(pat, ex)
    push!(items, result)

    items
end

function find_all_calls(ex::Expr)
    items = []
    no_arg_pat = :(_f())
    args_pat = :(_f(_x...))

    result = Espresso.findex(no_arg_pat, ex)
    push!(items, result...)

    result = Espresso.findex(args_pat, ex)
    push!(items, result...)

    for item in result
        if item isa Expr
            for arg in item.args[2:end]
                if arg isa Expr
                    push!(items, find_all_calls(arg)...)
                end
            end
        end
    end

    items
end

end

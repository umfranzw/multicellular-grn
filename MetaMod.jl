module MetaMod

function array_to_call(ex::Expr, f::Symbol)
    Expr(:call, f, ex.args...)
end

function call_to_array(ex::Expr)
    Expr(:vect, ex.args[2:end]...)
end

function array_to_foldl(ex::Expr, f::Symbol)
    Expr(:call, :foldl, f, ex)
end

function array_to_foldr(ex::Expr, f::Symbol)
    Expr(:call, :foldr, f, ex)
end

function fold_to_array(ex::Expr)
    ex.args[3]
end

function array_to_ref(ex::Expr, index::Int64)
    Expr(:ref, ex, index)
end

function ref_to_array(ex::Expr)
    ex.args[1]
end

function call_to_call(ex::Expr, f::Symbol)
    Expr(:call, f, ex.args[2:end]...)
end

function search(ex::Expr, key_ex::Expr)
    found = Array{Expr, 1}()
    if ex.head == key_ex.head
        args_present = true
        i = 1
        while args_present && i <= length(key_ex.args)
            key_arg = key_ex.args[i]
            args_present = key_arg in ex.args
            i += 1
        end

        if args_present
            push!(found, ex)
        end
    end

    for arg in ex.args
        if typeof(arg) == Expr
            push!(found, search(arg, key_ex)...)
        end
    end

    found
end

end

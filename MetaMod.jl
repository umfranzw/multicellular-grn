module MetaMod

import Espresso

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

end

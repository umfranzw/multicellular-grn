module SymMod

import Base.show

export Sym, SymType

@enum SymType::Int64 FcnCall FcnVal DataVar IntConst

struct Sym
    val::Any
    type::SymType
    arg_range::UnitRange #for fcns. Use -1 to indicate a variable number

    function Sym(val::Any, type::SymType, arg_range::UnitRange=0:0)
        new(val, type, arg_range)
    end

    function Sym(val::Any, type::SymType, arg_range::Int64)
        new(val, type, arg_range:arg_range)
    end
end

function show(io::IO, sym::Sym)
    print(io, sym.val)
end

end

module SymMod

export Sym, SymType

@enum SymType::Int64 FcnCall FcnVal DataVar IntConst

struct Sym
    val::Any
    type::SymType
end

end

module SymMod

export Sym, SymType

@enum SymType::Int64 FcnCall FcnVal DataVar IntConst

struct Sym
    val::Any
    type::SymType
    args_needed::Int64 #for fcns. Use -1 to indicate a variable number
end

end

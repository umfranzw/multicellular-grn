module SettingsMod

using SymMod

fcns = Set{Sym}()
terms = Set{Sym}()

#functions
push!(fcns,
      Sym(:+, SymMod.FcnCall, 1:-1),
      Sym(:-, SymMod.FcnCall, 1:-1),
      Sym(:*, SymMod.FcnCall, 1:-1),
      Sym(:÷, SymMod.FcnCall, 2),
      Sym(:^, SymMod.FcnCall, 2)
      )

push!(terms,
      Sym(:x, SymMod.DataVar),
      Sym(:1, SymMod.IntConst),
      Sym(:2, SymMod.IntConst),
      Sym(:3, SymMod.IntConst)
      )

end

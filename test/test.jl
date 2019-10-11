import RunMod
using IndividualMod
using SymMod
using CellMod
using CellTreeMod
using FitnessMod

run = RunMod.get_first_run()
#indiv = IndividualMod.rand_init(run)
#indiv.root_cell.sym = Sym(:*, SymMod.FcnCall, -1)

#three = Cell(run, indiv.genes, indiv.root_cell, Sym(:x, SymMod.DataVar, 0))
#mult = Cell(run, indiv.genes, indiv.root_cell, Sym(2, SymMod.IntConst, -1))

#FitnessMod.eval(indiv)

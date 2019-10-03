import RunMod
using CellMod
using GeneMod
using SymMod
using CellTreeMod

run = RunMod.get_first_run()
genes = map(i -> GeneMod.rand_init(run, i), 1:run.num_genes)
root = Cell(run, genes, nothing, Sym(:+, SymMod.FcnCall))

genes = map(i -> GeneMod.rand_init(run, i), 1:run.num_genes)
child1 = Cell(run, genes, root, nothing)

genes = map(i -> GeneMod.rand_init(run, i), 1:run.num_genes)
child2 = Cell(run, genes, root, Sym(2, SymMod.IntConst))

push!(root.children, child1, child2)

cells = CellTreeMod.find_empty(root)
println(cells)

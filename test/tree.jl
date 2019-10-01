import RunMod
using CellMod
using GeneMod

run = RunMod.get_first_run()
genes = map(i -> GeneMod.rand_init(run, i), 1:run.num_genes)
cell = Cell(run, genes, nothing, :x)

import RunMod
using IndividualMod
using GeneMod
using SymMod
using CellMod
using CellTreeMod
using ProteinStoreMod
using ProteinMod
using ProteinPropsMod

run = RunMod.get_first_run()
genes = map(i -> GeneMod.rand_init(run, i), 1:run.num_genes)

root = Cell(run, genes, nothing, Sym(:*, SymMod.FcnCall, -1))
left = Cell(run, genes, root, Sym(1, SymMod.IntConst, 0))
right = Cell(run, genes, root, Sym(2, SymMod.IntConst, 0))

p_root = Protein(
    run,
    ProteinProps(
        ProteinPropsMod.Reg,
        ProteinPropsMod.Intra,
        ProteinPropsMod.Activate,
        ProteinPropsMod.A
    ),
    zeros(run.num_genes)
)
ProteinStoreMod.insert(root.proteins, p_root, true)

p_left = Protein(
    run,
    ProteinProps(
        ProteinPropsMod.Reg,
        ProteinPropsMod.Intra,
        ProteinPropsMod.Activate,
        ProteinPropsMod.A
    ),
    zeros(run.num_genes)
)
ProteinStoreMod.insert(left.proteins, p_left, true)


p_right = Protein(
    run,
    ProteinProps(
        ProteinPropsMod.Reg,
        ProteinPropsMod.Intra,
        ProteinPropsMod.Activate,
        ProteinPropsMod.A
    ),
    zeros(run.num_genes)
)
ProteinStoreMod.insert(right.proteins, p_right, true)

import RunMod
using IndividualMod
using GeneMod
using SymMod
using CellMod
using CellTreeMod
using ProteinStoreMod
using ProteinMod
using ProteinPropsMod
using DiffusionMod

function test_diffusion()
    run = RunMod.get_first_run()
    genes = map(i -> GeneMod.rand_init(run, i), 1:run.num_genes)

    root = Cell(run, genes, nothing, Sym(:*, SymMod.FcnCall, -1))
    left = Cell(run, genes, root, Sym(1, SymMod.IntConst, 0))
    right = Cell(run, genes, root, Sym(2, SymMod.IntConst, 0))

    test_concs = zeros(run.num_genes)
    test_concs[run.num_genes ÷ 2] = 1.0

    p_root = Protein(
        run,
        ProteinProps(
            ProteinPropsMod.Reg,
            ProteinPropsMod.Intra,
            ProteinPropsMod.Activate,
            ProteinPropsMod.A
        ),
        copy(test_concs)
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
        copy(test_concs)
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
        copy(test_concs)
    )
    ProteinStoreMod.insert(right.proteins, p_right, true)

    p_inter = Protein(
        run,
        ProteinProps(
            ProteinPropsMod.Reg,
            ProteinPropsMod.Inter,
            ProteinPropsMod.Activate,
            ProteinPropsMod.A
        ),
        copy(test_concs)
    )
    ProteinStoreMod.insert(root.proteins, p_inter, true)

    #test inter-cellular diffusion
    DiffusionMod.diffuse_inter_cell_proteins(root)

    #make sure left child has the protein in an appropriate quantity
    inter_result = ProteinStoreMod.get(left.proteins, p_inter.props)
    @assert inter_result != nothing && (0 < maximum(inter_result.concs) < 1)

    #make sure the intra cell protein has not diffused
    intra_result = ProteinStoreMod.get(left.proteins, p_left.props)
    @assert intra_result != nothing && intra_result.concs == p_left.concs

    #test intra-cellular diffusion
    DiffusionMod.diffuse_intra_cell_proteins(root)
    #make sure the intra cell protein has diffused
    mid = run.num_genes ÷ 2
    intra_protein = ProteinStoreMod.get(root.proteins, p_root.props)
    @assert(intra_protein != nothing &&
            intra_protein.concs[mid] < test_concs[mid] &&
            intra_protein.concs[mid - 1] > 0 &&
            intra_protein.concs[mid + 1] > 0)
end

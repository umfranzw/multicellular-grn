module ProteinAppActionsMod

using ProteinPropsMod
using ProteinMod
using CellMod
using CellTreeMod
using SymMod
using GeneMod
using Printf

export AppAction

struct AppAction
    name::String
    fcn::Function
end

struct AppArgs
    tree::CellTree
    cell::Cell
    genes::Array{Gene, 1}
    protein::Protein
end

const app_actions = Array{AppAction, 1}([
    AppAction("make_parent_op_plus", args -> make_parent_op(args, :+)),
    AppAction("make_parent_op_minus", args -> make_parent_op(args, :-)),
    AppAction("make_parent_op_mult", args -> make_parent_op(args, :*)),
    AppAction("make_child_int_1", args -> make_child_int(args, 1)),
    AppAction("make_child_int_2", args -> make_child_int(args, 2)),
    AppAction("make_child_int_3", args -> make_child_int(args, 3)),
])

function push_app_actions()
    for action in app_actions
        push!(ProteinProps.AppActions, action)
    end
end

#returns a set containing any deleted cells
function make_parent_op(args::AppArgs, op::Symbol)
    if args.cell.parent == nothing
        #@info @sprintf("Applying protein: make_parent_op\n")
        parent = Cell(args.cell.config, args.genes, Sym(op, SymMod.FcnCall, -1))
        CellMod.add_parent(args.cell, parent)
        args.cell.energy /= 2
        args.tree.root = parent
    end

    Set{Cell}()
end

function make_child_int(args::AppArgs, val::Int64)
    if args.cell.sym.type == SymMod.FcnCall && (args.cell.sym.args_needed == -1 || length(args.cell.children) < args.cell.sym.args_needed)
        #@info @sprintf("Applying protein: make_child_int\n")
        child = Cell(args.cell.config, args.genes, Sym(val, SymMod.IntConst, 0))
        args.cell.energy /= 2
        CellMod.add_child(args.cell, child)
    end

    Set{Cell}()
end

function run_app_action(tree::CellTree, cell::Cell, genes::Array{Gene, 1}, protein::Protein)
    args = AppArgs(tree, cell, genes, protein)
    action = app_actions[Int64(protein.props.app_action)]
    action.fcn(args)
end

end

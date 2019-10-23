module ProteinAppActionsMod

using ProteinMod
using ProteinPropsMod
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

#returns a set containing any deleted cells
function make_parent_op(tree::CellTree, cell::Cell, genes::Array{Gene, 1}, protein::Protein)
    if cell.parent == nothing
        #@info @sprintf("Applying protein: make_parent_op\n")
        parent = Cell(cell.config, genes, Sym(:*, SymMod.FcnCall, -1))
        CellMod.add_parent(cell, parent)
        tree.root = parent
    end

    Set{Cell}()
end

function make_child_int(tree::CellTree, cell::Cell, genes::Array{Gene, 1}, protein::Protein)
    if cell.sym.type == SymMod.FcnCall && (cell.sym.args_needed == -1 || length(cell.children) < cell.sym.args_needed)
        #@info @sprintf("Applying protein: make_child_int\n")
        child = Cell(cell.config, genes, Sym(2, SymMod.IntConst, 0))
        CellMod.add_child(cell, child)
    end

    Set{Cell}()
end

const app_actions = Array{AppAction, 1}([
    AppAction("make_parent_op", make_parent_op),
    AppAction("make_child_int", make_child_int),
])

function run_app_action(tree::CellTree, cell::Cell, genes::Array{Gene, 1}, protein::Protein)
    action = app_actions[Int64(protein.props.app_action)]
    action.fcn(tree, cell, genes, protein)
end

end

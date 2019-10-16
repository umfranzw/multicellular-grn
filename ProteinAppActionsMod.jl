module ProteinAppActionsMod

using ProteinMod
using ProteinPropsMod
using CellMod
using SymMod

export AppAction

struct AppAction
    name::String
    fcn::Function
end

function make_parent_op(cell::Cell, protein::Protein)
    if cell.parent == nothing
        parent = Cell(cell.run, cell.genes, nothing, Sym(SymMod.FcnCall, :*, -1))
        push!(parent.children, cell)
    end
end

function make_child_int(cell::Cell, protein::Protein)
    if cell.sym.type == FcnCall && (cell.sym.args_needed == -1 || length(cell.children) < cell.sym.args_needed)
        child = Cell(cell.run, cell.genes, cell, Sym(SymMod.IntConst, 2, 0))
    end
end

const app_actions = Array{AppAction, 1}([
    AppAction("make_parent_op", make_parent_op),
    AppAction("make_child_int", make_child_int),
])

function run_app_action(cell::Cell, protein::Protein)
    action = app_actions[Int64(protein.props.app_action)]
    action.fcn(cell, protein)
end

end

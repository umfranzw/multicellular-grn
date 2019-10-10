module ProteinAppActionsMod

using ProteinMod
using ProteinPropsMod
using CellMod

export AppAction

struct AppAction
    #TODO: need to provide: daughter cell symbol, direction of app
    name::String
    fcn::Function
end

function temp(cell::Cell, protein::Protein)
end

const app_actions = Array{AppAction, 1}([
    AppAction("a", temp),
    AppAction("b", temp),
    AppAction("c", temp),
    AppAction("d", temp)
])

function get_app_action(protein::Protein)
    
end

function run_app_action(cell::Cell, protein::Protein)
    action = app_actions[Int64(protein.props.app_action)]
    action.fcn(cell, protein)
end

end

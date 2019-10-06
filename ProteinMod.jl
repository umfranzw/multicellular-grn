module ProteinMod

using RunMod
using BitUtilsMod
import RandUtilsMod
import Random
import MiscUtilsMod
import Base.copy
import ProteinAppActionsMod

export Protein, ProteinProps

@enum ProteinType::UInt8 Reg=0 App=1
@enum ProteinTarget::UInt8 Intra=0 Inter=1
@enum ProteinRegAction::UInt8 Activate=1 Inhibit=2
#note: these must match the values in the app_actions array
@enum ProteinAppAction::UInt8 A B C D

struct AppAction
    #TODO: need to provide: daughter cell symbol, direction of app
    name::String
    fcn::Function
end

const app_actions = Array{AppAction, 1}([
    AppAction("a", ProteinAppActionsMod.temp),
    AppAction("b", ProteinAppActionsMod.temp),
    AppAction("c", ProteinAppActionsMod.temp),
    AppAction("d", ProteinAppActionsMod.temp)
])

const num_types = MiscUtilsMod.num_enum_vals(ProteinType)
const num_targets = MiscUtilsMod.num_enum_vals(ProteinTarget)
const num_reg_actions = MiscUtilsMod.num_enum_vals(ProteinRegAction)
const num_app_actions = MiscUtilsMod.num_enum_vals(ProteinAppAction)

struct ProteinProps
    type::ProteinType
    target::ProteinTarget
    reg_action::ProteinRegAction
    app_action::ProteinAppAction #index of action in app_actions array
end

mutable struct Protein
    run::Run
    props::ProteinProps
    concs::Array{Float64, 1}

    function Protein(run::Run, props::ProteinProps, rand_concs::Bool)
        if rand_concs
            concs = RandUtilsMod.rand_floats(run, run.num_genes)
        else
            concs = zeros(Float64, run.num_genes)
        end
        
        new(run, props, concs)
    end

    function Protein(run::Run,  props::ProteinProps, concs::Array{Float64, 1})
        new(run, props, concs)
    end
end

function copy(protein::Protein)
    #only the concs need to be deep copied
    Protein(protein.run, protein.props, copy(protein.concs))
end

function get_app_action(protein::Protein)
    app_actions[Int64(protein.props.app_action)]
end

end

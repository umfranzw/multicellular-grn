module ProteinMod

using RunMod
using BitUtilsMod
import RandUtilsMod
import Random
import MiscUtilsMod
import Base.copy

export Protein

@enum ProteinType::UInt8 Reg=0 Growth=1 Inhibit=2
@enum ProteinTarget::UInt8 Intra=0 Inter=1
@enum ProteinRegAction::UInt8 RateUp=0 RateDown=1
#note: these must match the values in the growth_actions array
@enum ProteinGrowthAction::UInt8 A B C D

struct GrowthAction
    #TODO: need to provide: daughter cell symbol, direction of growth
    name::String
end

const growth_actions = Array{GrowthAction, 1}([
    GrowthAction("a"),
    GrowthAction("b"),
    GrowthAction("c"),
    GrowthAction("d")
])

const num_types = MiscUtilsMod.num_enum_vals(ProteinType)
const type_digits = MiscUtilsMod.digits_needed(num_types)
const num_targets = MiscUtilsMod.num_enum_vals(ProteinTarget)
const target_digits = MiscUtilsMod.digits_needed(num_targets)
const num_reg_actions = MiscUtilsMod.num_enum_vals(ProteinRegAction)
const reg_action_digits = MiscUtilsMod.digits_needed(num_reg_actions)
const num_growth_actions = MiscUtilsMod.num_enum_vals(ProteinGrowthAction)
const growth_action_digits = MiscUtilsMod.digits_needed(num_growth_actions)

mutable struct Protein
    run::Run
    type::ProteinType
    concs::Array{Float64, 1}
    target::ProteinTarget
    reg_action::ProteinRegAction
    growth_action::ProteinGrowthAction #index of action in growth_actions array

    function Protein(run::Run, rand_concs::Bool, type::ProteinType=Reg, target::ProteinTarget=Intra, reg_action::ProteinRegAction=RateUp, growth_action::ProteinGrowthAction=A)
        if rand_concs
            concs = RandUtilsMod.rand_floats(run, run.num_genes)
        else
            concs = zeros(Float64, run.num_genes)
        end
        
        new(run, type, concs, target, reg_action, growth_action)
    end

    function Protein(run::Run, concs::Array{Float64, 1}, type::ProteinType=Reg, target::ProteinTarget=Intra, reg_action::ProteinRegAction=RateUp, growth_action::ProteinGrowthAction=A)
        new(run, type, concs, target, reg_actions, growth_action)
    end
end

function copy(protein::Protein)
    #only the concs need to be deep copied
    Protein(protein.run, copy(protein.concs), protein.type, protein.target, protein.reg_action, protein.growth_action)
end

function get_growth_action(protein::Protein)
    get_growth_action(protein.growth_action)
end

function get_growth_action(action::ProteinGrowthAction)
    growth_actions[Int64(action)]
end

end

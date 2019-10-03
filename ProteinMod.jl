module ProteinMod

using RunMod
using BitUtilsMod
import RandUtilsMod
import Random
import MiscUtilsMod
import Base.copy

export Protein, ProteinProps

@enum ProteinType::UInt8 Reg=0 Growth=1
@enum ProteinTarget::UInt8 Intra=0 Inter=1
@enum ProteinRegAction::UInt8 Activate=1 Inhibit=2
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
# const type_digits = MiscUtilsMod.digits_needed(num_types)
# const type_fs = Formatting.FormatSpec("#0$(type_digits)d")

const num_targets = MiscUtilsMod.num_enum_vals(ProteinTarget)
# const target_digits = MiscUtilsMod.digits_needed(num_targets)
# const target_fs = Formatting.FormatSpec("#0$(target_digits)d")

const num_reg_actions = MiscUtilsMod.num_enum_vals(ProteinRegAction)
# const reg_action_digits = MiscUtilsMod.digits_needed(num_reg_actions)
# const reg_action_fs = Formatting.FormatSpec("#0$(reg_action_digits)d")

const num_growth_actions = MiscUtilsMod.num_enum_vals(ProteinGrowthAction)
# const growth_action_digits = MiscUtilsMod.digits_needed(num_growth_actions)
# const growth_action_fs = Formatting.FormatSpec("#0$(growth_action_digits)d")

struct ProteinProps
    type::ProteinType
    target::ProteinTarget
    reg_action::ProteinRegAction
    growth_action::ProteinGrowthAction #index of action in growth_actions array
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

function get_growth_action(protein::Protein)
    growth_actions[Int64(protein.props.growth_action)]
end

end

module ProteinPropsMod

import MiscUtilsMod
import Base.hash
import Base.==
import Base.copy

export ProteinProps,
    hash, ==

@enum ProteinType::UInt8 Reg=0 App=1
@enum ProteinTarget::UInt8 Intra=0 Inter=1
@enum ProteinRegAction::UInt8 Activate=0 Inhibit=1
#note: these must match the values in the app_actions array
@enum ProteinAppAction::UInt8 A B C D

const ProteinEnum = Union{ProteinType, ProteinTarget, ProteinRegAction, ProteinAppAction}

const num_types = MiscUtilsMod.num_enum_vals(ProteinType)
const num_targets = MiscUtilsMod.num_enum_vals(ProteinTarget)
const num_reg_actions = MiscUtilsMod.num_enum_vals(ProteinRegAction)
const num_app_actions = MiscUtilsMod.num_enum_vals(ProteinAppAction)

#note: must be only one field of each type
mutable struct ProteinProps
    type::ProteinType
    target::ProteinTarget
    reg_action::ProteinRegAction
    app_action::ProteinAppAction #index of action in app_actions array
end

function hash(props::ProteinProps)
    hash((props.type, props.target, props.reg_action, props.app_action))
end

==(p1::ProteinProps, p2::ProteinProps) = (p1.type == p2.type &&
                                          p1.target == p2.target &&
                                          p1.reg_action == p2.reg_action &&
                                          p1.app_action == p2.app_action)

function copy(props::ProteinProps)
    ProteinProps(props.type, props.target, props.reg_action, props.app_action)
end

end

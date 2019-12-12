module ProteinPropsMod

import Base.hash
import Base.==
import Base.copy
import Base.show
import Formatting
using MiscUtilsMod

export ProteinProps,
    hash, ==

@enum ProteinType::UInt8 Reg=0 App=1
@enum ProteinTarget::UInt8 Intra=0 Inter=1
@enum ProteinRegAction::UInt8 Activate=0 Inhibit=1
#note: these must match the values in the app_actions array in ProteinAppActionsMod
@enum ProteinAppAction::UInt8 P1=1 P2=2 P3=3 P4=4 P5=5 P6=6

const num_types = MiscUtilsMod.num_enum_vals(ProteinType)
const type_digits = MiscUtilsMod.digits_needed(num_types)
const type_fs = Formatting.FormatSpec("0$(type_digits)d")

const num_targets = MiscUtilsMod.num_enum_vals(ProteinTarget)
const target_digits = MiscUtilsMod.digits_needed(num_targets)
const target_fs = Formatting.FormatSpec("0$(target_digits)d")

const num_reg_actions = MiscUtilsMod.num_enum_vals(ProteinRegAction)
const reg_action_digits = MiscUtilsMod.digits_needed(num_reg_actions)
const reg_action_fs = Formatting.FormatSpec("0$(reg_action_digits)d")

const num_app_actions = MiscUtilsMod.num_enum_vals(ProteinAppAction)
const app_action_digits = MiscUtilsMod.digits_needed(num_app_actions)
const app_action_fs = Formatting.FormatSpec("0$(app_action_digits)d")

#note: must be only one field of each type
mutable struct ProteinProps
    type::ProteinType
    target::ProteinTarget
    reg_action::ProteinRegAction
    app_action::ProteinAppAction #index of action in app_actions array
end

function show(io::IO, props::ProteinProps, ilevel::Int64=0)
    pairs = (
        (type_fs, props.type),
        (target_fs, props.target),
        (reg_action_fs, props.reg_action),
        (app_action_fs, props.app_action)
    )
    for (fs, val) in pairs
        str = Formatting.fmt(fs, val)
        iprint(io, str, ilevel)
    end
    println(io, "")
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

function to_str(props::ProteinProps)
    buf = IOBuffer()
    show(buf, props)
    seek(buf, 0)

    chomp(read(buf, String)) #protein sequence string (remove the newline)
end

end

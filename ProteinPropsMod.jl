module ProteinPropsMod

import Base.hash
import Base.==
import Base.copy
import Base.show
import Formatting
import CustomEnumMod
using MiscUtilsMod

export ProteinProps,
    hash, ==

#@enum ProteinType::UInt8 Reg=0 App=1
CustomEnumMod.define_enum(
    :ProteinType,
    [:Reg, :App]
)
ProteinTypes = CustomEnumMod.ProteinTypes() #instance (contains enum values)
ProteinType = CustomEnumMod.ProteinTypeVal #value type (CustomEnumMod.ProteinTypeVal)

#@enum ProteinTarget::UInt8 Intra=0 Inter=1
CustomEnumMod.define_enum(
    :ProteinTarget,
    [:Intra, :Inter]
)
ProteinTargets = CustomEnumMod.ProteinTargets()
ProteinTarget = CustomEnumMod.ProteinTargetVal

#@enum ProteinRegAction::UInt8 Activate=0 Inhibit=1
CustomEnumMod.define_enum(
    :ProteinRegAction,
    [:Activate, :Inhibit]
)
ProteinRegActions = CustomEnumMod.ProteinRegActions()
ProteinRegAction = CustomEnumMod.ProteinRegActionVal

#note: these must match the values in the app_actions array in ProteinAppActionsMod
#@enum ProteinAppAction::UInt8 P1=1 P2=2 P3=3 P4=4 P5=5 P6=6
CustomEnumMod.define_enum(
    :ProteinAppAction,
    Symbol[] #values will be added dynamically by ProteinAppActionsMod
)
ProteinAppActions = CustomEnumMod.ProteinAppActions()
ProteinAppAction = CustomEnumMod.ProteinAppActionVal

#note: must be only one field of each type
mutable struct ProteinProps
    type::ProteinType
    target::ProteinTarget
    reg_action::ProteinRegAction
    app_action::ProteinAppAction
end

function show(io::IO, props::ProteinProps, ilevel::Int64=0)
    pairs = (
        (ProteinTypes, props.type),
        (ProteinTargets, props.target),
        (ProteinRegActions, props.reg_action),
        (ProteinAppActions, props.app_action)
    )
    for (enum, val) in pairs
        width = MiscUtilsMod.digits_needed(length(enum))
        fs = Formatting.FormatSpec("0$(width)d")
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

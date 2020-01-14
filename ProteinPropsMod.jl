module ProteinPropsMod

import Base.hash
import Base.==
#import Base.copy
import Base.show
import Formatting
using MiscUtilsMod

export ProteinProps,
    hash, ==

@enum ProteinType::Int8 Reg=1 App

@enum ProteinTarget::Int8 Intra=1 Inter

@enum ProteinRegAction::Int8 Activate=1 Inhibit

#note: these must match the values in the app_actions array in ProteinAppActionsMod
#@enum ProteinAppAction::UInt8 P1=1 P2=2 P3=3 P4=4 P5=5 P6=6
num_app_actions = 0

#note: must be only one field of each type
mutable struct ProteinProps
    type::ProteinType
    target::ProteinTarget
    reg_action::ProteinRegAction
    app_action::UInt8
end

function set_num_app_actions(num::Int64)
    global num_app_actions
    
    num_app_actions = UInt8(num)
end

function show(io::IO, props::ProteinProps, ilevel::Int64=0)
    global num_app_actions
    
    pairs = (
        (ProteinType, props.type),
        (ProteinTarget, props.target),
        (ProteinRegAction, props.reg_action)
    )
    for (enum, val) in pairs
        width = MiscUtilsMod.digits_needed(length(instances(enum)))
        fs = Formatting.FormatSpec("0$(width)d")
        str = Formatting.fmt(fs, Int64(val))
        iprint(io, str, ilevel)
    end
    #app_action
    width = MiscUtilsMod.digits_needed(Int64(num_app_actions))
    fs = Formatting.FormatSpec("0$(width)d")
    str = Formatting.fmt(fs, props.app_action)
    iprint(io, str, ilevel)
    
    println(io, "")
end

function hash(props::ProteinProps)
    hash((props.type, props.target, props.reg_action, props.app_action))
end

==(p1::ProteinProps, p2::ProteinProps) = (p1.type == p2.type &&
                                          p1.target == p2.target &&
                                          p1.reg_action == p2.reg_action &&
                                          p1.app_action == p2.app_action)

# function copy(props::ProteinProps)
#     ProteinProps(props.type, props.target, props.reg_action, props.app_action)
# end

function to_str(props::ProteinProps)
    buf = IOBuffer()
    show(buf, props)
    seek(buf, 0)

    chomp(read(buf, String)) #protein sequence string (remove the newline)
end

end

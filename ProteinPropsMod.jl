module ProteinPropsMod

import Base.hash
import Base.==
import Base.show
import Formatting
import Random
using RandUtilsMod
using MiscUtilsMod
using RunMod

export ProteinProps,
    hash, ==

@enum ProteinType::Int8 Internal=1 Neighbour Diffusion Application

@enum ProteinFcn::Int8 Inhibit=1 Activate

@enum ProteinAction::Int8 SymProb=1 Divide Sensor

@enum ProteinLoc::Int8 Top=1 Bottom Left Right

const opposite_locs = Dict{ProteinLoc, ProteinLoc}(
    Top => Bottom,
    Bottom => Top,
    Left => Right,
    Right => Left
)

#note: must be only one field of each type
mutable struct ProteinProps
    type::ProteinType
    fcn::ProteinFcn
    action::ProteinAction
    loc::ProteinLoc
    arg::Int8
end

function rand_init(
    config::Config;
    type::Union{Array{ProteinPropsMod.ProteinType, 1}, Nothing}=nothing,
    fcn::Union{Array{ProteinPropsMod.ProteinType, 1}, Nothing}=nothing,
    action::Union{Array{ProteinPropsMod.ProteinType, 1}, Nothing}=nothing,
    loc::Union{Array{ProteinPropsMod.ProteinLoc, 1}, Nothing}=nothing,
    arg::Union{UInt8, Nothing}=nothing
)
    pairs = (
        (type, ProteinType),
        (fcn, ProteinFcn),
        (action, ProteinAction),
        (loc, ProteinLoc)
    )

    vals = Array{Any, 1}()
    for (val, enum) in pairs
        if val == nothing
            instance = RandUtilsMod.rand_enum_val(config, enum)
        else
            instance = Random.rand(config.rng, val)
        end
        push!(vals, instance)
    end

    if arg == nothing
        arg_val = Random.rand(config.rng, Int8)
    else
        arg_val = arg
    end
    push!(vals, arg_val)

    ProteinProps(vals...)
end

function get_opposite_loc(loc::ProteinLoc)
    global opposite_locs
    
    opposite_locs[loc]
end

function show(io::IO, props::ProteinProps, ilevel::Int64=0)
    pairs = (
        (ProteinType, props.type),
        (ProteinFcn, props.fcn),
        (ProteinAction, props.action),
        (ProteinLoc, props.loc)
    )
    for (enum, val) in pairs
        # width = MiscUtilsMod.digits_needed(length(instances(enum)))
        # fs = Formatting.FormatSpec("0$(width)d")
        # str = Formatting.fmt(fs, Int64(val))
        # iprint(io, str, ilevel)
        iprint(io, string(val)[1:3], ilevel) #print first 3 chars of value name
    end
    iprint(io, string(props.arg), ilevel)
    
    println(io, "")
end

function hash(props::ProteinProps)
    hash((props.type, props.fcn, props.action, props.loc, props.arg))
end

==(p1::ProteinProps, p2::ProteinProps) = (p1.type == p2.type &&
                                          p1.fcn == p2.fcn &&
                                          p1.loc == p2.loc &&
                                          p1.action == p2.action &&
                                          p1.arg == p2.arg)

function to_str(props::ProteinProps)
    buf = IOBuffer()
    show(buf, props)
    seek(buf, 0)

    chomp(read(buf, String)) #protein sequence string (remove the newline)
end

end

module ProteinPropsMod

import Base.hash
import Base.==
import Base.show
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
    arg::UInt8
end

function rand_prop(
    config::Config,
    enum::Type{T},
    vals::Union{Array{T, 1}, Nothing}=nothing
) where {T <: Union{ProteinPropsMod.ProteinType,
                    ProteinPropsMod.ProteinFcn,
                    ProteinPropsMod.ProteinAction,
                    ProteinPropsMod.ProteinLoc,
                    UInt8}}
    if vals == nothing
        if enum == UInt8
            instance = Random.rand(config.rng, UInt8)
        else
            instance = RandUtilsMod.rand_enum_val(config, enum)
        end
    else
        instance = Random.rand(config.rng, vals)
    end

    instance
end

function rand_init(
    config::Config;
    type::Union{Array{ProteinPropsMod.ProteinType, 1}, Nothing}=nothing,
    fcn::Union{Array{ProteinPropsMod.ProteinFcn, 1}, Nothing}=nothing,
    action::Union{Array{ProteinPropsMod.ProteinAction, 1}, Nothing}=nothing,
    loc::Union{Array{ProteinPropsMod.ProteinLoc, 1}, Nothing}=nothing,
    arg::Union{Array{UInt8, 1}, Nothing}=nothing
)
    pairs = (
        (ProteinType, type),
        (ProteinFcn, fcn),
        (ProteinAction, action),
        (ProteinLoc, loc),
        (UInt8, arg)
    )
    vals = map(p -> rand_prop(config, p...), pairs)
    
    ProteinProps(vals...)
end

function get_opposite_loc(loc::ProteinLoc)
    global opposite_locs
    
    opposite_locs[loc]
end

function get_abbrev_props_str(;
                              type::Union{ProteinPropsMod.ProteinType, Nothing}=nothing,
                              fcn::Union{ProteinPropsMod.ProteinFcn, Nothing}=nothing,
                              action::Union{ProteinPropsMod.ProteinAction, Nothing}=nothing,
                              loc::Union{ProteinPropsMod.ProteinLoc, Nothing}=nothing
                              )
    buf = IOBuffer()
    for val in (type, fcn, action, loc)
        if val != nothing
            chunk = string(val)[1:3] #print first 3 chars of value name
            print(buf, chunk)
        end
    end

    String(take!(buf))
end

function show(io::IO, props::ProteinProps, ilevel::Int64=0)
    desc = get_abbrev_props_str(
        type=props.type,
        fcn=props.fcn,
        action=props.action,
        loc=props.loc
    )

    iprint(io, desc, ilevel)
    iprint(io, "-", ilevel)
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

function to_str(props::ProteinProps, is_initial::Bool)
    buf = IOBuffer()
    show(buf, props)
    if is_initial
        seek(buf, buf.size - sizeof("\n"))
        write(buf, "*\n")
    end
    seek(buf, 0)

    string(chomp(read(buf, String))) #protein sequence string (remove the newline)
end

end

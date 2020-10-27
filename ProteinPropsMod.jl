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

@enum ProteinType::Int8 Internal Neighbour Diffusion Application

@enum ProteinFcn::Int8 Inhibit=-1 Activate=1

@enum ProteinAction::Int8 SymProb Divide #Sensor

@enum ProteinLoc::Int8 Right Top Left #Bottom...

#note: must be only one field of each type
mutable struct ProteinProps
    type::ProteinType
    tag::UInt8
    action::ProteinAction
    arg::Int8
end

function get_fcn(props::ProteinProps)
    get_fcn(props.arg)
end

function get_fcn(arg::Int8)
    arg < 0 ? Inhibit : Activate
end

function rand_prop(
    config::Config,
    enum::Type{T},
    vals::Union{Array{T, 1}, Nothing}=nothing
) where {T <: Union{ProteinPropsMod.ProteinType, ProteinPropsMod.ProteinAction, Int8, UInt8}}
    if vals == nothing
        #choose a completely random option
        if enum == UInt8
            instance = Random.rand(config.rng, UInt8(0):UInt8(config.run.tag_limit))
        elseif enum == Int8
            instance = Random.rand(config.rng, enum)
        else
            instance = RandUtilsMod.rand_enum_val(config, enum)
        end
    else
        #choose from one of the provided vals
        instance = Random.rand(config.rng, vals)
    end

    instance
end

function rand_init(
    config::Config;
    type::Union{Array{ProteinPropsMod.ProteinType, 1}, Nothing}=nothing,
    tag::Union{Array{UInt8, 1}, Nothing}=nothing,
    action::Union{Array{ProteinPropsMod.ProteinAction, 1}, Nothing}=nothing,
    fcn::Union{ProteinPropsMod.ProteinFcn, Nothing}=nothing,
    arg::Union{Array{ProteinPropsMod.ProteinAction, 1}, Array{ProteinPropsMod.ProteinLoc, 1}, Nothing}=nothing
)
    type_val = rand_prop(config, ProteinType, type)
    tag_val = rand_prop(config, UInt8, tag)
    action_val = rand_prop(config, ProteinAction, action)
    arg_val = rand_prop(config, Int8, arg)
    if fcn != nothing
        #force sign to be fcn
        arg_val = Int8(fcn) * abs(arg_val)
    end
    
    ProteinProps(type_val, tag_val, action_val, arg_val)
end

function get_abbrev_props_str(;
                              type::Union{ProteinPropsMod.ProteinType, Nothing}=nothing,
                              tag::Union{UInt8, Nothing}=nothing,
                              action::Union{ProteinPropsMod.ProteinAction, Nothing}=nothing,
                              arg::Union{Int8, Nothing}=nothing
                              )
    labels = Array{String, 1}()
    if type != nothing
        push!(labels, string(type)[1:3])
    end
    if tag != nothing
        push!(labels, string(tag))
    end
    if action != nothing
        push!(labels, string(action)[1:3])
    end
    if arg != nothing
        push!(labels, string(arg))
    end
    
    join(labels, ":")
end

function show(io::IO, props::ProteinProps, ilevel::Int64=0)
    desc = get_abbrev_props_str(
        type=props.type,
        tag=props.tag,
        action=props.action,
        arg=props.arg
    )

    iprint(io, desc, ilevel)
    println(io, "")
end

function hash(props::ProteinProps)
    hash((props.type, props.tag, props.action, props.arg))
end

==(p1::ProteinProps, p2::ProteinProps) = (p1.type == p2.type &&
                                          p1.tag == p2.tag &&
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

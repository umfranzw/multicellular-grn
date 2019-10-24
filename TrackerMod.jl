module TrackerMod

import CodecZlib
using RunMod
using IndividualMod
using Serialization
using CellTreeMod
using Printf

@enum BestType RunBest GenBest

export Tracker

tracker = nothing
run_best = nothing
gen_best = nothing

mutable struct Tracker
    run::Run
    ea_states::Dict{String, Array{Array{UInt8, 1}, 1}}
    reg_states::Dict{String, Array{Array{UInt8, 1}, 1}}

    function Tracker(run::Run)
        new(run, Dict{String, Array{Array{UInt8, 1}, 1}}(), Dict{String, Array{Array{UInt8, 1}, 1}}())
    end
end

function create_tracker(run::Run)
    global tracker
    tracker = Tracker(run)
end

function get_tracker()
    global tracker
    tracker
end

function destroy_tracker()
    global tracker
    tracker = nothing
end

function update_bests(pop::Array{Individual, 1})
    global tracker
    global gen_best
    global run_best
    
    rb_updated = false
    gb_updated = false
    for indiv in pop
        if gen_best == nothing || indiv.fitness < gen_best.fitness
            gen_best = indiv
            gb_updated = true
            
            if run_best == nothing || indiv.fitness < run_best.fitness
                run_best = indiv
                rb_updated = true
                
            end
        end
    end

    if gb_updated
        @info join(
            (
                "gen_best:",
                @sprintf("fitness: %0.2f", gen_best.fitness),
                CellTreeMod.to_expr_str(gen_best.cell_tree)
            ),
            "\n"
        )
        
        if rb_updated
            @info join(
                (
                    "run_best:",
                    @sprintf("fitness: %0.2f", run_best.fitness),
                    CellTreeMod.to_expr_str(run_best.cell_tree)
                ),
                "\n"
            )
        end
    end
end

function save_ea_state(label::String, pop::Array{T, 1}) where T
    global tracker
    
    buf = IOBuffer()
    Serialization.serialize(buf, pop)
    #we'll compress the data *inside* the dictionary. This makes it small enough that it's reasonable to keep it
    #in memory (avoiding a disk access) during the simulation
    bytes = CodecZlib.transcode(CodecZlib.GzipCompressor, buf.data)
    
    if label ∉ keys(tracker.ea_states)
        tracker.ea_states[label] = Array{Array{UInt8, 1}, 1}()
    end
    push!(tracker.ea_states[label], bytes)
end

function save_ea_state_on_step(ea_iter::Int64, label::String, pop::Array{T, 1}) where T
    global tracker

    if ea_iter in tracker.run.step_range
        save_ea_state(label, pop)
    end
end

function save_reg_state(label::String, pop_trees::Array{Array{CellTree, 1}, 1})
    global tracker

    buf = IOBuffer()
    Serialization.serialize(buf, pop_trees)
    bytes = CodecZlib.transcode(CodecZlib.GzipCompressor, buf.data)

    if label ∉ keys(tracker.reg_states)
        tracker.reg_states[label] = Array{Array{UInt8, 1}, 1}()
    end
    push!(tracker.reg_states[label], bytes)
end

function write_data(path::String)
    global tracker

    #note: it's really not worth compressing the enclosing dictionaries here - since the data inside is already compressed,
    #it's small enough...
    out_stream = open(path, "w")
    Serialization.serialize(out_stream, (tracker.ea_states, tracker.reg_states))
    close(out_stream)
end

end

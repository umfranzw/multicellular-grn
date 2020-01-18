module TrackerMod

import CodecZlib
using RunMod
using IndividualMod
using Serialization
using CellTreeMod
using Printf

@enum BestType::UInt8 RunBest GenBest
@enum TagType::UInt8 EAState RegState RunState

export Tracker

tracker = nothing

#in bytes (max that CodecZlib.GzipCompressor can handle is 2^32 - 1)
#note: 2^20 bytes = 1MB
const cache_size = 512 * 2^20

mutable struct Tracker
    run::Run
    path::String
    file_handle::IOStream
    run_best::Union{Individual, Nothing}
    gen_best::Union{Individual, Nothing}
    cache::IOBuffer
end

function create_tracker(run::Run, path::String)
    global tracker

    file_handle = open(path, "w")
    cache = IOBuffer(maxsize=cache_size)
    
    tracker = Tracker(run, path, file_handle, nothing, nothing, cache)
    save_run()

    tracker
end

function get_tracker()
    global tracker
    tracker
end

function destroy_tracker()
    global tracker

    flush_cache()
    close(tracker.file_handle)
    tracker = nothing
end

function update_bests(pop::Array{Individual, 1})
    global tracker
    
    rb_updated = false
    gb_updated = false
    for indiv in pop
        if tracker.gen_best == nothing || indiv.fitness < tracker.gen_best.fitness
            tracker.gen_best = deepcopy(indiv)
            gb_updated = true
            
            if tracker.run_best == nothing || indiv.fitness < tracker.run_best.fitness
                tracker.run_best = tracker.gen_best #can just use the copy that was already made
                rb_updated = true
                
            end
        end
    end

    if gb_updated
        # @info join(
        #     (
        #         "gen_best:",
        #         @sprintf("fitness: %0.2f", tracker.gen_best.fitness),
        #         CellTreeMod.to_expr_str(tracker.gen_best.cell_tree)
        #     ),
        #     "\n"
        # )
        
        if rb_updated
            @info join(
                (
                    "run_best:",
                    @sprintf("fitness: %0.2f", tracker.run_best.fitness),
                    CellTreeMod.to_expr_str(tracker.run_best.cell_tree)
                ),
                "\n"
            )
        end
    end
end

function write_to_cache(data::Array{UInt8, 1})
    global tracker

    bytes_needed = sizeof(Int64) + length(data) #we will need to write the size (Int64) and then data
    if tracker.cache.size > cache_size
        flush_cache()
    end
    write(tracker.cache, Int64(length(data)))
    write(tracker.cache, data)
end

function flush_cache()
    global tracker

    println("flushing cache")
    bytes = CodecZlib.transcode(CodecZlib.GzipCompressor, take!(tracker.cache))

    #write size, then compressed data
    write(tracker.file_handle, Int64(length(bytes)))
    write(tracker.file_handle, bytes)
end

function save_run()
    global tracker

    if tracker.run.log_data
        state = (RunState, tracker.run)
        buf = IOBuffer()
        Serialization.serialize(buf, state)
        write_to_cache(take!(buf))
    end
end

function save_ea_state(pop::Array{Individual, 1}, ea_step::Int64, force::Bool=false)
    for i in 1:length(pop)
        save_ea_state(pop[i], ea_step, i, force)
    end
end

function save_ea_state(indiv::Individual, ea_step::Int64, index::Int64, force::Bool=false)
    global tracker

    if tracker.run.log_data && (ea_step in tracker.run.step_range || force)
        state = (EAState, ea_step, index, indiv)
        buf = IOBuffer()
        Serialization.serialize(buf, state)
        write_to_cache(take!(buf))
    end
end

function save_reg_state(tree::CellTree, ea_step::Int64, reg_step::Int64, index::Int64)
    global tracker

    if tracker.run.log_data
        state = (RegState, ea_step, reg_step, index, tree)
        buf = IOBuffer()
        Serialization.serialize(buf, state)
        write_to_cache(take!(buf))
    end
end

end

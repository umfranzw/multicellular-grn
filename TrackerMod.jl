module TrackerMod

import CodecXz
using RunMod
using IndividualMod
using Serialization
using CellTreeMod
using BestInfoMod
using CheckpointMod
using Printf

@enum StateType::UInt8 IndivState RunState RunBestInfoState FitnessesState
@enum StateTime::UInt8 AfterBind AfterProd
@enum StateContentType::UInt8 Changes Checkpoint

export Tracker

tracker = nothing

#in bytes (max that CodecZlib.GzipCompressor can handle is 2^32 - 1)
#note: 2^20 bytes = 1MB
#const cache_size = 512 * 2^20

checkpoint_reg_steps = Set{Int64}([1])

mutable struct LastCheckpoint
    indiv_bytes::Array{UInt8, 1}
    ea_step::Int64
    pop_index::Int64
    reg_step::Int64
end

mutable struct Tracker
    run::Run
    path::String
    file_handle::Union{IOStream, Nothing}
    run_best::BestInfo
    gen_best::BestInfo
    fitnesses::Array{Array{Float64, 1}, 1}
    file_handle_lock::ReentrantLock
    last_checkpoint::Union{LastCheckpoint, Nothing}
end

function create_tracker(run::Run, path::String)
    global tracker

    if run.log_level > RunMod.LogNone
        file_handle = open(path, "w")
    else
        file_handle = nothing
    end
    tracker = Tracker(run, path, file_handle, BestInfo(), BestInfo(), Array{Array{Float64, 1}, 1}(), ReentrantLock(), nothing)
    save_run()

    tracker
end

function get_tracker()
    global tracker
    tracker
end

function destroy_tracker()
    global tracker

    if tracker.run.log_level > RunMod.LogNone
        close(tracker.file_handle)
    end
    tracker = nothing
end

function update_fitnesses(pop::Array{Individual, 1}, ea_step::Int64)
    global tracker
    
    rb_updated = false
    gb_updated = false
    fitnesses = Array{Float64, 1}()
    for pop_index in 1:length(pop)
        indiv = pop[pop_index]
        push!(fitnesses, indiv.fitness)
        
        if update(tracker.gen_best, indiv, ea_step, pop_index, tracker.run.reg_steps + 1)
            gb_updated = true
            #note: we only ever need to update run best if gen best was updated
            if update(tracker.run_best, tracker.gen_best.indiv, ea_step, pop_index, tracker.run.reg_steps + 1, make_copy=false) #can just use the copy that was already made for gen best
                rb_updated = true
            end
        end
    end

    if rb_updated
        @info join(
            (
                "run_best:",
                @sprintf("fitness: %0.2f", tracker.run_best.indiv.fitness),
                CellTreeMod.to_expr_str(tracker.run_best.indiv.cell_tree)
            ),
            "\n"
        )
    end
    push!(tracker.fitnesses, fitnesses)
end

function write_obj(state_type::StateType, step_tag::Array{Int64, 1}, obj_bytes::Array{UInt8, 1}, state_time::Union{StateTime, Nothing}=nothing, state_content_type::Union{StateContentType, Nothing}=nothing)
    global tracker

    buf = IOBuffer()
    compressed = CodecXz.transcode(CodecXz.XzCompressor, take!(buf))
    
    lock(tracker.file_handle_lock)
    write(tracker.file_handle, UInt8(state_tag)) #always write state type
    #optionally write the other two
    if state_time != nothing
        write(tracker.file_handle, UInt8(state_time))
    end
    if state_content_type != nothing
        write(tracker.file_handle, UInt8(state_content_type))
    end

    #always write step_tag (ea_step, indiv_index, reg_step)
    for tag_val in step_tag
        write(tracker.file_handle, tag_val)
    end
    write(tracker.file_handle, Int64(length(compressed))) #write size
    write(tracker.file_handle, compressed) #write obj
    unlock(tracker.file_handle_lock)
end

function save_run()
    global tracker

    if tracker.run.log_level >= RunMod.LogFitnesses
        buf = IOBuffer()
        Serialization.serialize(buf, tracker.run)
        write_obj(RunState, Array{Int64, 1}(), take!(run))
    end
end

function save_run_best()
    global tracker

    if tracker.run.log_level >= RunMod.LogIndivs
        buf = IOBuffer()
        Serialization.serialize(buf, tracker.run_best)
        write_obj(RunBestInfoState, Array{Int64, 1}(), take!(buf))
    end
end

function save_fitnesses()
    global tracker

    if tracker.run.log_level >= RunMod.LogFitnesses
        buf = IOBuffer()
        Serialization.serialize(buf, tracker.fitnesses)
        write_obj(FitnessesState, Array{Int64, 1}(), take!(buf))
    end
end

function save_reg_state(indiv::Individual, ea_step::Int64, reg_step::Int64, pop_index::Int64, state_time::StateTime)
    global tracker

    if tracker.run.log_level >= RunMod.LogIndivs
        if reg_step ∈ TrackerMod.checkpoint_reg_steps
            buf = IOBuffer()
            Serialization.serialize(buf, indiv)
            data = take!(buf)
            state_content_type = Checkpoint
            tracker.last_checkpoint = LastCheckpoint(data, ea_step, pop_index, reg_step)
        else
            if tracker.last_checkpoint == nothing
                throw Exception("Tracker error: Attempted to save reg state before any checkpoints have been stored")
            elseif tracker.last_checkpoint.ea_step == ea_step && tracker.last_checkpoint.pop_index == pop_index && tracker.last_checkpoint.reg_step != reg_step
                buf = IOBuffer()
                Serialization.serialize(buf, indiv)
                change_info = CheckpointMod.compress(tracker.last_checkpoint.indiv_bytes, take!(buf))
                Serialization.serialize(buf, change_info)
                data = take!(buf)
                state_content_type = Changes
            else
                throw Exception("Tracker error: Attempted to save state using invalid last checkpoint")
            end
        end
        
        
        write_obj(IndivState, Array{Int64, 1}([ea_step, reg_step, pop_index]), data, state_time, state_content_type)
    end
end

end

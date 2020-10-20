module TrackerMod

using RunMod
using IndividualMod
using Serialization
using CellTreeMod
using BestInfoMod
using CheckpointMod
using Printf
using DataStructures
using CompressionMod
import Statistics

@enum StateType::UInt8 IndivState RunState RunBestInfoState FitnessesState
@enum StateTime::UInt8 AfterBind AfterProd
@enum StateContentType::UInt8 Changes Checkpoint

export Tracker

tracker = nothing

#in bytes (max that CodecZlib.GzipCompressor can handle is 2^32 - 1)
#note: 2^20 bytes = 1MB
#const cache_size = 512 * 2^20

#Note: items should be inserted in ascending order
checkpoint_reg_steps = OrderedSet{Int64}([1])
file_chunk_size = 1024 #bytes

mutable struct LastCheckpoint
    indiv_bytes::Array{UInt8, 1}
    ea_step::Int64
    pop_index::Int64
    reg_step::Int64
end

mutable struct Tracker
    run::Run
    path::String
    local_files::Union{Array{IOStream, 1}, Nothing}
    run_best::BestInfo
    gen_best::BestInfo
    fitnesses::Array{Array{Float64, 1}, 1}
    last_checkpoint::Union{Array{Union{LastCheckpoint, Nothing}, 1}, Nothing}
end

function create_tracker(run::Run, path::String)
    global tracker

    if run.log_level > RunMod.LogNone
        local_files = Array{IOStream, 1}()
        last_checkpoint = Array{Union{LastCheckpoint, Nothing}, 1}()
        for i in 1:Threads.nthreads()
            local_fname = get_local_fname(path, i)
            push!(local_files, open(local_fname, "w+")) #open for writing *and reading* (see destroy_tracker())
            push!(last_checkpoint, nothing)
        end
    else
        local_files = nothing
        last_checkpoint = nothing
    end
    
    tracker = Tracker(run, path, local_files, BestInfo(), BestInfo(), Array{Array{Float64, 1}, 1}(), last_checkpoint)
    save_compression_type()
    save_run()

    tracker
end

function get_tracker()
    global tracker
    tracker
end

function get_local_fname(path::String, threadid::Int64)
    "$(path)-$(threadid)"
end

function destroy_tracker()
    global tracker

    size = -1
    if tracker.run.log_level > RunMod.LogNone
        #copy each of the tracker.local_files into a single global file (location is given by tracker.path)
        #note: it's important that local_files[1] is written first, since it contains the compression_type and run, and they need to be at the top (DataMod expects it)
        global_file = open(tracker.path, "w")
        for i in 1:Threads.nthreads()
            info = stat(tracker.local_files[i])
            seek(tracker.local_files[i], 0)
            while !eof(tracker.local_files[i])
                chunk = read(tracker.local_files[i], TrackerMod.file_chunk_size)
                write(global_file, chunk)
            end
            close(tracker.local_files[i])
            fname = get_local_fname(tracker.path, i)
            rm(fname)
        end
        close(global_file)
        file_info = stat(tracker.path)
        size = file_info.size #in bytes
    end
    tracker = nothing

    size
end

#note: this will always be called with in a single threaded context, so there's no need for locking
function update_fitnesses(pop::Array{Individual, 1}, ea_step::Int64, output_buf::IOBuffer)
    global tracker
    
    rb_updated = false
    gb_updated = false
    fitnesses = Array{Float64, 1}()
    for pop_index in 1:length(pop)
        indiv = pop[pop_index]
        push!(fitnesses, indiv.fitness)
        
        if BestInfoMod.update(tracker.gen_best, indiv, ea_step, pop_index, tracker.run.reg_steps + 1)
            gb_updated = true
            #note: we only ever need to update run best if gen best was updated
            if BestInfoMod.update(tracker.run_best, tracker.gen_best.indiv, ea_step, pop_index, tracker.run.reg_steps + 1, make_copy=false) #can just use the copy that was already made for gen best
                rb_updated = true
            end
        end
    end

    mean = Statistics.mean(fitnesses)
    std = Statistics.std(fitnesses, mean=mean)
    write(output_buf, @sprintf("mean fitness: %0.5f\n", mean))
    write(output_buf, @sprintf("std fitness: %0.5f\n", std))
    
    if rb_updated
        write(output_buf, "run_best updated:\n")
        write(output_buf, @sprintf("\tfitness: %0.5f\n", tracker.run_best.indiv.fitness))
        write(output_buf, @sprintf("\texpr: %s\n", CellTreeMod.to_expr_str(tracker.run_best.indiv.cell_tree)))
    end
    push!(tracker.fitnesses, fitnesses)
end

function write_obj(
    state_type::StateType,
    comp_obj_bytes::Array{UInt8, 1},
    step_tag::Union{Array{Int64, 1}, Nothing}=nothing,
    state_time::Union{StateTime, Nothing}=nothing,
    state_content_type::Union{StateContentType, Nothing}=nothing
)
    global tracker

    #format: state_type, size, compressed_data, [step_tag], [state_time], [state_content_type]
    
    #this needs to be thread safe
    #lock(tracker.file_handle_lock)

    file_index = Threads.threadid()
    
    #always write state type, size, and compressed data
    #println("writing state_type: $(state_type)")
    write(tracker.local_files[file_index], UInt8(state_type))
    #println("writing size: $(length(comp_obj_bytes))")
    write(tracker.local_files[file_index], Int64(length(comp_obj_bytes)))
    #println("writing comp_obj_bytes")
    write(tracker.local_files[file_index], comp_obj_bytes)
    
    #optionally write the others
    #step tag format: (ea_step, indiv_index, reg_step)
    if step_tag != nothing
        #println("writing step tag: $(step_tag)")
        for tag_val in step_tag
            write(tracker.local_files[file_index], tag_val)
        end
    end

    if state_time != nothing
        #println("writing state_time: $(state_time)")
        write(tracker.local_files[file_index], UInt8(state_time))
    end

    if state_content_type != nothing
        #println("writing state_content_type: $(state_content_type)")
        write(tracker.local_files[file_index], UInt8(state_content_type))
    end
    
    #unlock(tracker.file_handle_lock)
end

function save_compression_type()
    global tracker

    if tracker.run.log_level > RunMod.LogNone
        write(tracker.local_files[Threads.threadid()], UInt8(tracker.run.compression_alg))
    end
end

function save_run()
    global tracker

    if tracker.run.log_level > RunMod.LogNone
        #println("saving run")
        buf = IOBuffer()
        Serialization.serialize(buf, tracker.run)
        compressed = CompressionMod.compress(tracker.run.compression_alg, take!(buf))
        write_obj(RunState, compressed)
        # println()
    end
end

function save_run_best()
    global tracker

    if tracker.run.log_level >= RunMod.LogIndivs
        #println("saving run_best")
        buf = IOBuffer()
        Serialization.serialize(buf, tracker.run_best)
        compressed = CompressionMod.compress(tracker.run.compression_alg, take!(buf))
        write_obj(RunBestInfoState, compressed)
        #println()
    end
end

function save_fitnesses()
    global tracker

    if tracker.run.log_level >= RunMod.LogFitnesses
        #println("saving fitnesses")
        buf = IOBuffer()
        Serialization.serialize(buf, tracker.fitnesses)
        compressed = CompressionMod.compress(tracker.run.compression_alg, take!(buf))
        write_obj(FitnessesState, compressed)
        #println()
    end
end

#This function should be thread-safe as long as each thread sticks to it's indiv, doing all reg_steps sequentially
function save_reg_state(indiv::Individual, ea_step::Int64, pop_index::Int64, reg_step::Int64, state_time::StateTime)
    global tracker

    #println("Saving reg state for: ($(ea_step), $(pop_index), $(reg_step))")
    
    if tracker.run.log_level >= RunMod.LogIndivs
        file_index = Threads.threadid()
        #note: only AfterBind indivs can be checkpoints - all AfterProd indivs use the last AfterBind indiv as their checkpoint base for checkpoint compression
        if state_time == AfterBind && reg_step ∈ TrackerMod.checkpoint_reg_steps
            #println("It's a checkpoint indiv")
            buf = IOBuffer()
            Serialization.serialize(buf, indiv)
            data = CompressionMod.compress(tracker.run.compression_alg, take!(buf))
            state_content_type = Checkpoint

            tracker.last_checkpoint[file_index] = LastCheckpoint(data, ea_step, pop_index, reg_step)
        else
            if tracker.last_checkpoint == nothing
                throw(Exception("Tracker error: Attempted to save reg state before any checkpoints have been stored"))
            elseif tracker.last_checkpoint[file_index].ea_step == ea_step && tracker.last_checkpoint[file_index].pop_index == pop_index && (tracker.last_checkpoint[file_index].reg_step < reg_step || state_time == AfterProd) #note: all AfterProd indivs use the last AfterBind checkpoint indiv as their checkpoint base for checkpoint compression
                #println("It's a change indiv")
                buf = IOBuffer()
                Serialization.serialize(buf, indiv)
                compressed = CompressionMod.compress(tracker.run.compression_alg, take!(buf))
                change_info = CheckpointMod.compress(tracker.last_checkpoint[file_index].indiv_bytes, compressed)
                Serialization.serialize(buf, change_info)
                data = CompressionMod.compress(tracker.run.compression_alg, take!(buf))
                state_content_type = Changes
            else
                throw(ErrorException("Tracker error: Attempted to save state using invalid last checkpoint"))
            end
        end

        
        # buf = IOBuffer()
        # Serialization.serialize(buf, indiv)
        # data = CodecXz.transcode(CodecXz.XzCompressor, take!(buf))
        # state_content_type = Checkpoint
        
        write_obj(IndivState, data, Array{Int64, 1}([ea_step, pop_index, reg_step]), state_time, state_content_type)

        #println()
    end
end

end

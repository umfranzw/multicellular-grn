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
    save_compression_type()
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
        
        if BestInfoMod.update(tracker.gen_best, indiv, ea_step, pop_index, tracker.run.reg_steps + 1)
            gb_updated = true
            #note: we only ever need to update run best if gen best was updated
            if BestInfoMod.update(tracker.run_best, tracker.gen_best.indiv, ea_step, pop_index, tracker.run.reg_steps + 1, make_copy=false) #can just use the copy that was already made for gen best
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
    lock(tracker.file_handle_lock)

    #always write state type, size, and compressed data
    #println("writing state_type: $(state_type)")
    write(tracker.file_handle, UInt8(state_type))
    #println("writing size: $(length(comp_obj_bytes))")
    write(tracker.file_handle, Int64(length(comp_obj_bytes)))
    #println("writing comp_obj_bytes")
    write(tracker.file_handle, comp_obj_bytes)
    
    #optionally write the others
    #step tag format: (ea_step, indiv_index, reg_step)
    if step_tag != nothing
        #println("writing step tag: $(step_tag)")
        for tag_val in step_tag
            write(tracker.file_handle, tag_val)
        end
    end

    if state_time != nothing
        #println("writing state_time: $(state_time)")
        write(tracker.file_handle, UInt8(state_time))
    end

    if state_content_type != nothing
        #println("writing state_content_type: $(state_content_type)")
        write(tracker.file_handle, UInt8(state_content_type))
    end
    
    unlock(tracker.file_handle_lock)
end

function save_compression_type()
    global tracker

    write(tracker.file_handle, UInt8(tracker.run.compression_alg))
end

function save_run()
    global tracker

    if tracker.run.log_level >= RunMod.LogFitnesses
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

function save_reg_state(indiv::Individual, ea_step::Int64, pop_index::Int64, reg_step::Int64, state_time::StateTime)
    global tracker

    #println("Saving reg state for: ($(ea_step), $(pop_index), $(reg_step))")

    if tracker.run.log_level >= RunMod.LogIndivs
        #note: only AfterBind indivs can be checkpoints - all AfterProd indivs use the last AfterBind indiv as their checkpoint base for checkpoint compression
        if state_time == AfterBind && reg_step ∈ TrackerMod.checkpoint_reg_steps
            #println("It's a checkpoint indiv")
            buf = IOBuffer()
            Serialization.serialize(buf, indiv)
            data = CompressionMod.compress(tracker.run.compression_alg, take!(buf))
            state_content_type = Checkpoint
            
            tracker.last_checkpoint = LastCheckpoint(data, ea_step, pop_index, reg_step)
        else
            if tracker.last_checkpoint == nothing
                throw(Exception("Tracker error: Attempted to save reg state before any checkpoints have been stored"))
            elseif tracker.last_checkpoint.ea_step == ea_step && tracker.last_checkpoint.pop_index == pop_index && (tracker.last_checkpoint.reg_step < reg_step || state_time == AfterProd) #note: all AfterProd indivs use the last AfterBind checkpoint indiv as their checkpoint base for checkpoint compression
                #println("It's a change indiv")
                buf = IOBuffer()
                Serialization.serialize(buf, indiv)
                compressed = CompressionMod.compress(tracker.run.compression_alg, take!(buf))
                change_info = CheckpointMod.compress(tracker.last_checkpoint.indiv_bytes, compressed)
                Serialization.serialize(buf, change_info)
                data = CompressionMod.compress(tracker.run.compression_alg, take!(buf))
                state_content_type = Changes
            else
                throw(Exception("Tracker error: Attempted to save state using invalid last checkpoint"))
            end
        end

        
        # buf = IOBuffer()
        # Serialization.serialize(buf, indiv)
        # data = CodecXz.transcode(CodecXz.XzCompressor, take!(buf))
        # state_content_type = Checkpoint
        
        write_obj(IndivState, data, Array{Int64, 1}([ea_step, pop_index, reg_step]), state_time, state_content_type)

        if ea_step == 1 && pop_index == 6 && reg_step ∈ (1, 2)
            creating = !isfile("/home/wayne/Documents/school/thesis/multicellular-grn/data/test")
            test_handle = open("/home/wayne/Documents/school/thesis/multicellular-grn/data/test", "a")
            temp = tracker.file_handle
            tracker.file_handle = test_handle
            if creating
                save_compression_type()
                save_run()
            end
            write_obj(IndivState, data, Array{Int64, 1}([ea_step, pop_index, reg_step]), state_time, state_content_type)
            tracker.file_handle = temp
            close(test_handle)
        end
            
        #println()
    end
end

end

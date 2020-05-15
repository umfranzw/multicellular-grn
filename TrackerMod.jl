module TrackerMod

import CodecZlib
using RunMod
using IndividualMod
using Serialization
using CellTreeMod
using Printf

@enum BestType::UInt8 RunBest GenBest
@enum TagType::UInt8 IndivState RunState RunBestInfoState FitnessesState

export Tracker, BestInfo

tracker = nothing

#in bytes (max that CodecZlib.GzipCompressor can handle is 2^32 - 1)
#note: 2^20 bytes = 1MB
#const cache_size = 512 * 2^20

mutable struct BestInfo
    index::Union{Nothing, Tuple{Int64, Int64, Int64}} #(ea_step, pop_index, reg_step) - reg step will always be run.reg_steps + 1
    indiv::Union{Nothing, Individual}

    function BestInfo()
        new(nothing, nothing)
    end
end

function update(info::BestInfo, indiv::Individual, ea_step::Int64, pop_index::Int64, reg_step::Int64; make_copy::Bool=true)
    updated = false
    if !is_set(info) || indiv.fitness < info.indiv.fitness
        if make_copy
            info.indiv = deepcopy(indiv)
        else
            info.indiv = indiv
        end
        
        info.index = (ea_step, pop_index, reg_step)
        updated = true
    end

    updated
end

function is_set(info::BestInfo)
    info.indiv != nothing
end

mutable struct Tracker
    run::Run
    path::String
    file_handle::Union{IOStream, Nothing}
    run_best::BestInfo
    gen_best::BestInfo
    fitnesses::Array{Array{Float64, 1}, 1}
    file_handle_lock::ReentrantLock
end

function create_tracker(run::Run, path::String)
    global tracker

    if run.log_level > RunMod.LogNone
        file_handle = open(path, "w")
    else
        file_handle = nothing
    end
    tracker = Tracker(run, path, file_handle, BestInfo(), BestInfo(), Array{Array{Float64, 1}, 1}(), ReentrantLock())
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

function write_obj(tag_type::TagType, tag::Array{Int64, 1}, obj::Any)
    global tracker

    buf = IOBuffer()
    Serialization.serialize(buf, obj)
    compressed = CodecZlib.transcode(CodecZlib.GzipCompressor, buf.data)
    
    lock(tracker.file_handle_lock)
    write(tracker.file_handle, UInt8(tag_type)) #write type
    for tag_val in tag
        write(tracker.file_handle, tag_val) #write tag
    end
    write(tracker.file_handle, Int64(length(compressed))) #write size
    write(tracker.file_handle, compressed) #write obj
    unlock(tracker.file_handle_lock)
end

function save_run()
    global tracker

    if tracker.run.log_level >= RunMod.LogFitnesses
        write_obj(RunState, Array{Int64, 1}(), tracker.run)
    end
end

function save_run_best()
    global tracker

    if tracker.run.log_level >= RunMod.LogIndivs
        write_obj(RunBestInfoState, Array{Int64, 1}(), tracker.run_best)
    end
end

function save_fitnesses()
    global tracker

    if tracker.run.log_level >= RunMod.LogFitnesses
        write_obj(FitnessesState, Array{Int64, 1}(), tracker.fitnesses)
    end
end

# function save_ea_state(pop::Array{Individual, 1}, ea_step::Int64, force::Bool=false)
#     global tracker
    
#     if tracker.run.log_data
#         for i in 1:length(pop)
#             if ea_step in tracker.run.step_range || force
#                 write_obj(IndivState, Array{Int64, 1}([ea_step, i]), pop[i])
#             end
#         end
#     end
# end

function save_reg_state(indiv::Individual, ea_step::Int64, reg_step::Int64, index::Int64)
    global tracker

    if tracker.run.log_level >= RunMod.LogIndivs
        write_obj(IndivState, Array{Int64, 1}([ea_step, reg_step, index]), indiv)
    end
end

end

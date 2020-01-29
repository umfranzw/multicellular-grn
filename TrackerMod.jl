module TrackerMod

import CodecZlib
using RunMod
using IndividualMod
using Serialization
using CellTreeMod
using Printf

@enum BestType::UInt8 RunBest GenBest
@enum TagType::UInt8 IndivState CellTreeState RunState

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
end

function create_tracker(run::Run, path::String)
    global tracker

    file_handle = open(path, "w")
    tracker = Tracker(run, path, file_handle, nothing, nothing)
    save_run()

    tracker
end

function get_tracker()
    global tracker
    tracker
end

function destroy_tracker()
    global tracker

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

function write_obj(tag_type::TagType, tag::Array{Int64, 1}, obj::Any)
    global tracker
    
    buf = IOBuffer()
    Serialization.serialize(buf, obj)
    compressed = CodecZlib.transcode(CodecZlib.GzipCompressor, buf.data)
    write(tracker.file_handle, UInt8(tag_type)) #write type
    for tag_val in tag
        write(tracker.file_handle, tag_val) #write tag
    end
    write(tracker.file_handle, Int64(length(compressed))) #write size
    write(tracker.file_handle, compressed) #write obj
end

function save_run()
    global tracker

    if tracker.run.log_data
        write_obj(RunState, Array{Int64, 1}(), tracker.run)
    end
end

function save_ea_state(pop::Array{Individual, 1}, ea_step::Int64, force::Bool=false)
    global tracker
    
    if tracker.run.log_data
        for i in 1:length(pop)
            if ea_step in tracker.run.step_range || force
                write_obj(IndivState, Array{Int64, 1}([ea_step, i]), pop[i])
            end
        end
    end
end

function save_reg_state(tree::CellTree, ea_step::Int64, reg_step::Int64, index::Int64)
    global tracker

    if tracker.run.log_data
        write_obj(CellTreeState, Array{Int64, 1}([ea_step, reg_step, index]), tree)
    end
end

end

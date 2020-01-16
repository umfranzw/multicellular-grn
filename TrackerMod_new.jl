module TrackerMod

import CodecZlib
using RunMod
using IndividualMod
using Serialization
using CellTreeMod
using Printf

@enum BestType::UInt8 RunBest GenBest
@enum TagType::UInt8 EAState RegState

export Tracker

tracker = nothing
run_best = nothing
gen_best = nothing

mutable struct Tracker
    run::Run
    path::String
    file_handle::IOStream
end

function create_tracker(run::Run, path::String)
    global tracker

    file_handle = open(path, "w")
    tracker = Tracker(run, path)
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
    global gen_best
    global run_best
    
    rb_updated = false
    gb_updated = false
    for indiv in pop
        if gen_best == nothing || indiv.fitness < gen_best.fitness
            gen_best = deepcopy(indiv)
            gb_updated = true
            
            if run_best == nothing || indiv.fitness < run_best.fitness
                run_best = gen_best #can just use the copy that was already made
                rb_updated = true
                
            end
        end
    end

    if gb_updated
        # @info join(
        #     (
        #         "gen_best:",
        #         @sprintf("fitness: %0.2f", gen_best.fitness),
        #         CellTreeMod.to_expr_str(gen_best.cell_tree)
        #     ),
        #     "\n"
        # )
        
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

function save_run()
    global tracker

    if tracker.run.log_data
        run_buf = IOBuffer()
        Serialization.serialize(run_buf, tracker.run)

        write(tracker.file_handle, Int64(run_buf.size))
        write(tracker.file_handle, run_buf.data)
    end
end

function save_ea_state(indiv::Individual, ea_step::Int64, index::Int64)
    global tracker

    if tracker.run.log_data
        indiv_buf = IOBuffer()
        Serialization.serialize(indiv_buf, indiv)
        
        #write tag
        write(tracker.file_handle, UInt8(EAState))
        write(tracker.file_handle, Int64(ea_step))
        write(tracker.file_handle, Int64(0)) #reg_step
        write(tracker.file_handle, Int64(index))
        write(tracker.file_handle, Int64(indiv__buf.size))

        #write individual
        write(tracker.file_handle, indiv_buf.data)
    end
end

function save_reg_state(tree::CellTree, ea_step::Int64, reg_step::Int64, index::Int64)
    global tracker

    if tracker.run.log_data
        tree_buf = IOBuffer()
        Serialization.serialize(tree_buf, tree)

        #write tag
        write(tracker.file_handle, UInt8(RegState))
        write(tracker.file_handle, Int64(ea_step))
        write(tracker.file_handle, Int64(reg_step))
        write(tracker.file_handle, Int64(index))
        write(tracker.file_handle, Int64(tree_buf.size))

        #write tree
        write(tracker.file_handle, tree_buf.data)
    end
end

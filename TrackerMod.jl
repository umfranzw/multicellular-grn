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
    states::Array{Array{UInt8, 1}}

    function Tracker(run::Run)
        new(run, Array{Tuple{Int64, Array{UInt8, 1}}, 1}())
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

function track_state(iter::Int64, pop::Array{Individual, 1})
    global tracker
    global gen_best
    global run_best
    
    if iter in tracker.run.step_range
        for indiv in pop
            if gen_best == nothing || indiv.fitness < gen_best.fitness
                gen_best = indiv
                info_str = join(
                    (
                        "gen_best:",
                        @sprintf("fitness: %0.2f", gen_best.fitness),
                        CellTreeMod.to_expr_str(gen_best.cell_tree)
                    ),
                    "\n"
                )
                @info info_str
                
                if run_best == nothing || indiv.fitness < run_best.fitness
                    run_best = indiv
                    info_str = join(
                        (
                            "run_best:",
                            @sprintf("fitness: %0.2f", run_best.fitness),
                            CellTreeMod.to_expr_str(run_best.cell_tree)
                        ),
                        "\n"
                    )
                    @info info_str
                end
            end
        end

        entry = (iter, pop)
        buf = IOBuffer()
        Serialization.serialize(buf, entry)
        bytes = CodecZlib.transcode(CodecZlib.GzipCompressor, buf.data)
        push!(tracker.states, bytes)
    end
end

function write_data(path::String)
    global tracker
    
    out_stream = open(path, "w")
    
    for state in tracker.states
        len = length(state)
        out_stream.write(Int64(len))
        out_stream.write(state)
    end
    
    close(out_stream)
end

end

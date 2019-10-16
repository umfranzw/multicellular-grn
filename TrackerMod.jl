module TrackerMod

using RunMod
using IndividualMod

@enum BestType RunBest GenBest

export Tracker

tracker = nothing

mutable struct Tracker
    run::Run
    states::Array{Tuple{Int64, Array{Individual, 1}}} #pairs of the form (iter, pop)

    function Tracker(run::Run)
        tracker = new(run, Array{Tuple{Int64, Array{Individual, 1}}, 1}())
    end
end


function create_tracker(run::Run)
    tracker = Tracker(run)
end

function get_tracker()
    tracker
end

function destroy_tracker()
    tracker = nothing
end

function track_state(iter::Int64, pop::Array{Individual, 1})
end

end

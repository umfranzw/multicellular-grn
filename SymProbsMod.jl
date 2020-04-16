module SymProbsMod

using SymMod
using SettingsMod
using RunMod
using RandUtilsMod
using MiscUtilsMod
using Printf
import Base.length
import Base.show

export SymProbs

index_to_sym = vcat([x for x in values(SettingsMod.fcns)], [x for x in values(SettingsMod.terms)])
sym_to_index = Dict{Sym, Int64}(zip(index_to_sym, 1:length(index_to_sym)))

mutable struct SymProbs
    probs::Dict{Sym, Float64}

    function SymProbs()
        initial_probs = repeat([1.0 / length(index_to_sym)], length(index_to_sym))
        new(Dict{Sym, Float64}(zip(index_to_sym, initial_probs)))
    end
end

function length(probs::SymProbs)
    length(probs.probs)
end

function alter_prob(probs::SymProbs, index::Int64, amount::Float64)
    global index_to_sym
    
    sym = index_to_sym[index]
    cur_prob = probs.probs[sym]
    new_prob = clamp(cur_prob + amount, 0.0, 1.0)
    probs.probs[sym] = new_prob

    #normalize so the sum stays equal to one
    total = sum(probs.probs)
    probs.probs ./= total
end

function choose_sym(probs::SymProbs, config::Config)
    #normalize the probabilities and build the "roulette wheel"
    total = foldl((s, p) -> s + p, values(probs.probs); init=0.0)
    next_bound = 0.0
    wheel = Array{Array{Union{Sym, Float64}, 1}, 1}() #[(Sym, cummulative prob), ...]
    for (sym, prob) in probs.probs
        next_bound += prob / total
        push!(wheel, [sym, next_bound])
    end

    wheel[end][2] = 1.0 #just in case we've got some floating point error

    sel_index = 1
    r = RandUtilsMod.rand_float(config) #random value in [0, 1)
    while r >= wheel[sel_index][2]
        sel_index += 1
    end

    wheel[sel_index][1] #return the selected symbol
end

function show(io::IO, probs::SymProbs, ilevel::Int64=0)
    #sort from highest to lowest
    pairs = [(sym, prob) for (sym, prob) in probs.probs]
    sort!(pairs, by=p -> p[2], rev=true)

    iprintln(io, "SymProbs:", ilevel)
    for (sym, prob) in pairs
        prob = @sprintf("%0.2f", prob)
        iprintln(io, "$(sym): $(prob)", ilevel + 1)
    end
end

end

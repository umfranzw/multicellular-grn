module SymProbsMod

using SymMod
using SettingsMod
using RunMod
import Base.length

export SymProbs

index_to_sym = vcat([x for x in values(SettingsMod.fcns)], [x for x in values(SettingsMod.terms)])
sym_to_index = Dict{Sym, Int64}(zip(index_to_sym, 1:length(index_to_sym)))

mutable struct SymProbs
    probs::Dict{Sym, Float64}

    function SymProbs()
        new(Dict{Sym, Float64}(zip(index_to_sym, zeros(length(index_to_sym)))))
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
end

function choose_sym(probs::SymProbs, config::Config)
    #normalize the probabilities and build the "roulette wheel"
    total = foldl((s, p) -> s + p, values(probs.probs); init=0.0)
    next_bound = 0.0
    wheel = Array{Tuple{Sym, Float64}}() #[(Sym, cummulative prob), ...]
    for (sym, prob) in probs.probs
        next_bound += prob / total
        push!(wheel, (sym, next_bound))
    end

    wheel[end] = 1.0 #just in case we've got some floating point error

    sel_index = 1
    r = RandUtilsMod.rand_float(gs.config) #random value in [0, 1)
    while r >= wheel[sel_index][2]
        sel_index += 1
    end

    wheel[sel_index][1] #return the selected symbol
end

end

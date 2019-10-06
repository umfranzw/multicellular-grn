module FitnessMod

using IndividualMod

function eval(indiv::Individual)
    indiv.fitness = 1.0
end

end

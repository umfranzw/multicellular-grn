module EvAlgMod

using RunMod
using IndividualMod
using MutateMod

function ev_alg(run::Run)
    pop = init_pop(run)

    for ea_step in 1:run.ea_steps
        MutateMod.mutate(run, pop)
        RegSimMod.reg_sim(run, pop)
        #TODO: reset cell trees and bindings before next ea_step
        #also need to be able to re-init cell (since initial proteins may have changed)
    end
end

function init_pop(run::Run)
    map(i -> IndividualMod.rand_init(run), 1:run.pop_size)
end

using RunMod
using EvAlgMod

function main()
    for run in RunMod.get_run_channel()
        EvAlgMod.ev_alg(run)
    end
end

main()

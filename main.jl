using RunMod
using EvAlgMod
using ProteinAppActionsMod

function main()
    ProteinAppActionsMod.push_app_actions()
    
    for run in RunMod.get_run_channel()
        EvAlgMod.ev_alg(run)
    end
end

main()
